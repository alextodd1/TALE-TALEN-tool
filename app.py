from flask import Flask, render_template, request, jsonify
from flask_wtf import FlaskForm
from wtforms import HiddenField, SelectField, TextAreaField, SubmitField, IntegerField
from wtforms.validators import DataRequired, Regexp, NumberRange, Length, ValidationError, Optional
import re
import time
from database_setup import db, TALEPair, init_db, generate_short_id


app = Flask(__name__)
app.secret_key = 'your_secret_key'
init_db(app)

talen_code_for_RVD = {"C": "HD", "T": "NG", "A": "NI", "G": "NN"}

class DNASequenceForm(FlaskForm):
    dna_sequence = TextAreaField('DNA Sequence', validators=[
        DataRequired(),
        Regexp(r'^[ATCGatcg]*$', message="Invalid DNA sequence. Only A, T, C, and G are allowed."),
        Length(max=100_000, message="DNA sequence must be 100,000 bases or less.") 
    ])
    min_tale_length = IntegerField('Minimum TALE Length', validators=[
        NumberRange(min=10, max=30, message="Minimum TALE length must be between 10 and 30.")
    ], default=18)
    max_tale_length = IntegerField('Maximum TALE Length', validators=[
        NumberRange(min=10, max=30, message="Maximum TALE length must be between 10 and 30.")
    ], default=18)
    min_spacer_length = IntegerField('Minimum Spacer Length', validators=[
        NumberRange(min=1, max=100, message="Minimum spacer length must be between 1 and 100.")
    ], default=25)
    max_spacer_length = IntegerField('Maximum Spacer Length', validators=[
        NumberRange(min=1, max=100, message="Maximum spacer length must be between 1 and 100.")
    ], default=30)
    search_position = IntegerField('Search Position', validators=[
        Optional(),
        NumberRange(min=1, message="Search position must be a positive integer.")
    ])
    search_range = IntegerField('Search Range (±)', validators=[
        Optional(),
        NumberRange(min=1, message="Search range must be a positive integer.")
    ])


    talen_code_g = SelectField('TALE Code for G', choices=[('NN', 'NN'), ('NH', 'NH')], default='NH')
    submit = SubmitField('Submit DNA Sequence')
    hidden_talen_code_g = HiddenField('Hidden TALE Code for G')

    def validate_max_spacer_length(form, field):
        if field.data - form.min_spacer_length.data > 30:
            raise ValidationError('The difference between maximum and minimum spacer length must not exceed 30.')
        
    def validate_max_tale_length(form, field):
        if field.data - form.min_tale_length.data > 10:
            raise ValidationError('The difference between maximum and minimum TALE length must not exceed 10.')


def generate_complementary_dna(sequence):
    complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return ''.join(complement_map[base] for base in sequence)

def calculate_cg_content(sequence):
    cg_count = sum(1 for base in sequence if base in 'CG')
    return (cg_count / len(sequence)) * 100

def is_cpg_island(sequence, start, end):
    length = end - start
    if length < 200:
        return False
    
    gc_count = sequence[start:end].count('G') + sequence[start:end].count('C')
    gc_percentage = (gc_count / length) * 100
    
    cpg_count = sequence[start:end].count('CG')
    expected_cpg = (sequence[start:end].count('C') * sequence[start:end].count('G')) / length
    observed_to_expected_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
    
    return gc_percentage > 50 and observed_to_expected_ratio > 0.6

def group_consecutive_bases(bases):
    ranges = []
    start = end = None
    for base in sorted(bases):
        if start is None:
            start = end = base
        elif base == end + 1:
            end = base
        else:
            ranges.append((start, end))
            start = end = base
    if start is not None:
        ranges.append((start, end))
    return ranges


def find_tale_pairs(sequence, min_tale_length=18, max_tale_length=18, min_spacer_length=25, max_spacer_length=30, search_position=None, search_range=None):
    comp_seq = generate_complementary_dna(sequence)
    pairs = []
    consecutive_at_regex = re.compile(r'[AT]{7,}')
    discarded_count = 0
    gc_discarded_count = 0
    strong_rvd_discarded_count = 0
    cpg_island_discarded_count = 0
    cpg_island_bases = set()

    if search_position is not None and search_range is not None:
        start_boundary = max(1, search_position - search_range)
        end_boundary = min(len(sequence), search_position + search_range)
    else:
        start_boundary = 1
        end_boundary = len(sequence)

    # Pre-compute GC content for all possible tale regions
    gc_content = [0] * (len(sequence) + 1)
    for i in range(1, len(sequence) + 1):
        gc_content[i] = gc_content[i-1] + (sequence[i-1] in 'GC')

    def calculate_gc_percentage(start, end):
        return ((gc_content[end] - gc_content[start]) / (end - start)) * 100

    def count_strong_rvds(rvd_sequence):
        return sum(1 for rvd in zip(rvd_sequence[::2], rvd_sequence[1::2]) if rvd in [('N', 'N'), ('H', 'D')])

    def is_cpg_island(start, end):
        if end - start < 200:
            return False
        gc_percentage = calculate_gc_percentage(start, end)
        cpg_count = sequence[start:end].count('CG')
        expected_cpg = ((gc_content[end] - gc_content[start]) / 2) ** 2 / (end - start)
        observed_to_expected_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
        return gc_percentage > 50 and observed_to_expected_ratio > 0.6

    for tale_length in range(min_tale_length, max_tale_length + 1):
        for i in range(start_boundary - 1, end_boundary - tale_length):
            if sequence[i] == 'T':
                tale_region = sequence[i+1 : i+1+tale_length]
                
                if is_cpg_island(max(0, i-100), min(len(sequence), i+tale_length+100)):
                    cpg_island_discarded_count += 1
                    cpg_island_bases.update(range(max(0, i-100), min(len(sequence), i+tale_length+100)))

                    continue

                if consecutive_at_regex.search(tale_region):
                    discarded_count += 1
                    continue

                if calculate_gc_percentage(i+1, i+1+tale_length) < 25:
                    gc_discarded_count += 1
                    continue

                rvd_sequence = ''.join(talen_code_for_RVD[base] for base in tale_region)
                start_pos = i
                end_pos = i + tale_length

                for spacer_length in range(min_spacer_length, min(max_spacer_length, 100) + 1):
                    comp_start = end_pos + spacer_length + 1
                    comp_end = comp_start + tale_length

                    if comp_start >= 0 and comp_end <= len(comp_seq) - 1 and comp_seq[comp_end] == 'T':
                        comp_tale_region = comp_seq[comp_start : comp_end]

                        if is_cpg_island(max(0, comp_start-100), min(len(comp_seq), comp_end+100)):
                            cpg_island_discarded_count += 1
                            continue

                        if consecutive_at_regex.search(comp_tale_region):
                            discarded_count += 1
                            continue

                        if calculate_gc_percentage(comp_start, comp_end) < 25:
                            gc_discarded_count += 1
                            continue

                        comp_rvd_sequence = ''.join(talen_code_for_RVD[base] for base in comp_tale_region[::-1])

                        if count_strong_rvds(rvd_sequence) < 3 or count_strong_rvds(comp_rvd_sequence) < 3:
                            strong_rvd_discarded_count += 1
                            continue

                        pairs.append((start_pos + 1, end_pos + 1, comp_start - 1, comp_end - 1, rvd_sequence, comp_rvd_sequence, spacer_length, tale_length))

    pairs.sort(key=lambda x: x[0])
    cpg_island_ranges = group_consecutive_bases(cpg_island_bases)


    return pairs, discarded_count, gc_discarded_count, strong_rvd_discarded_count, cpg_island_discarded_count, cpg_island_ranges


@app.route('/', methods=['GET', 'POST'])
def index():
    form = DNASequenceForm()
    if request.method == 'POST' and form.validate_on_submit():
        start_time = time.time()

        dna_sequence = form.dna_sequence.data.upper()
        min_tale_length = form.min_tale_length.data
        max_tale_length = form.max_tale_length.data
        min_spacer_length = form.min_spacer_length.data
        max_spacer_length = form.max_spacer_length.data
        g_code = form.talen_code_g.data
        search_position = form.search_position.data
        search_range = form.search_range.data
        talen_code_for_RVD["G"] = g_code
        form.hidden_talen_code_g.data = form.talen_code_g.data 

        if search_position and search_range:
            start_boundary = max(1, search_position - search_range)
            end_boundary = min(len(dna_sequence), search_position + search_range)
            searched_range = f"{start_boundary}-{end_boundary}"
        else:
            searched_range = None        
                
        comp_dna_sequence = generate_complementary_dna(dna_sequence)
        tale_pairs, discarded_count, gc_discarded_count, strong_rvd_discarded_count, cpg_island_discarded_count, cpg_island_ranges = find_tale_pairs(dna_sequence, min_tale_length, max_tale_length, min_spacer_length, max_spacer_length, search_position, search_range)         
        session_id = generate_short_id()

        # Prepare bulk insert data
        tale_pair_objects = [
            TALEPair(
                session_id=session_id,
                start=pair[0],
                end=pair[1],
                comp_start=pair[2],
                comp_end=pair[3],
                rvd=pair[4],
                comp_rvd=pair[5],
                spacer_length=pair[6],
                tale_length=pair[7],
                g_code=g_code
            ) for pair in tale_pairs
        ]

        # Perform bulk insert
        db.session.bulk_save_objects(tale_pair_objects)
        db.session.commit()

        tale_pairs_count = len(tale_pairs)
        execution_time = time.time() - start_time

        return render_template('index.html', form=form, dna_sequence=dna_sequence, comp_dna_sequence=comp_dna_sequence, 
                               min_tale_length=min_tale_length, max_tale_length=max_tale_length, 
                               min_spacer_length=min_spacer_length, max_spacer_length=max_spacer_length, 
                               discarded_count=discarded_count, gc_discarded_count=gc_discarded_count,
                               strong_rvd_discarded_count=strong_rvd_discarded_count,
                               cpg_island_discarded_count=cpg_island_discarded_count, 
                               cpg_island_ranges=cpg_island_ranges,
                               dna_sequence_length=len(dna_sequence), tale_pairs_count=tale_pairs_count,
                               execution_time=execution_time, session_id=session_id, searched_range=searched_range)
    return render_template('index.html', form=form)

@app.route('/api/tale_pairs', methods=['GET'])
def get_tale_pairs():
    session_id = request.args.get('session_id')
    page = request.args.get('page', 1, type=int)
    per_page = 100

    tale_pairs = TALEPair.query.filter_by(session_id=session_id).paginate(page=page, per_page=per_page, error_out=False)
    
    return jsonify({
        'tale_pairs': [
            {
                'start': pair.start,
                'end': pair.end,
                'comp_start': pair.comp_start,
                'comp_end': pair.comp_end,
                'rvd': pair.rvd,
                'comp_rvd': pair.comp_rvd,
                'spacer_length': pair.spacer_length,
                'tale_length': pair.tale_length,
                'g_code': pair.g_code
            } for pair in tale_pairs.items
        ],
        'total': tale_pairs.total,
        'pages': tale_pairs.pages,
        'current_page': tale_pairs.page
    })


@app.route('/about')
def about():
    return render_template('about.html')

if __name__ == '__main__':
    app.run(debug=False)

    
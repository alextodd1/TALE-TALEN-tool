from flask_sqlalchemy import SQLAlchemy
from datetime import datetime, timezone
import string
import random

db = SQLAlchemy()

def generate_short_id(length=6):
    """Generate a short random string of specified length."""
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

class TALEPair(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    session_id = db.Column(db.String(6), index=True)
    start = db.Column(db.Integer)
    end = db.Column(db.Integer)
    comp_start = db.Column(db.Integer)
    comp_end = db.Column(db.Integer)
    rvd = db.Column(db.String(60))
    comp_rvd = db.Column(db.String(60))
    spacer_length = db.Column(db.Integer)
    tale_length = db.Column(db.Integer)
    g_code = db.Column(db.String(2))
    created_at = db.Column(db.DateTime, default=datetime.now(timezone.utc))

def init_db(app):
    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///tale_pairs.db'
    app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
    db.init_app(app)
    
    with app.app_context():
        db.create_all()
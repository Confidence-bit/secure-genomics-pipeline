from cryptography.fernet import Fernet
from logging_module import log

def load_key():
    with open("fernet.key", "rb") as f:
        return f.read()

def encrypt_file(input_file, output_file):
    key = load_key()
    cipher = Fernet(key)

    with open(input_file, "rb") as f:
        data = f.read()

    encrypted = cipher.encrypt(data)

    with open(output_file, "wb") as f:
        f.write(encrypted)

    log(f"Encrypted: {input_file} → {output_file}")

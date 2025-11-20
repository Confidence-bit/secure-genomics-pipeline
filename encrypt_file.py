from cryptography.fernet import Fernet
import sys

key = open("fernet.key", "rb").read()
cipher = Fernet(key)

with open(sys.argv[1], "rb") as f:
    plaintext = f.read()

ciphertext = cipher.encrypt(plaintext)

with open(sys.argv[2], "wb") as f:
    f.write(ciphertext)

print("File encrypted successfully!")

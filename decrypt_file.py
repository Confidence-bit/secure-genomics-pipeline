from cryptography.fernet import Fernet
import sys

key = open("fernet.key", "rb").read()
cipher = Fernet(key)

with open(sys.argv[1], "rb") as f:
    ciphertext = f.read()

plaintext = cipher.decrypt(ciphertext)

with open(sys.argv[2], "wb") as f:
    f.write(plaintext)

print("File decrypted successfully!")

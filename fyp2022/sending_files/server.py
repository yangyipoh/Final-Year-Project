import socket

IP = "192.168.1.102"
PORT = 12345
ADDR = (IP, PORT)
SIZE = 1024
FORMAT = "utf-8"

def main():
	print("[STARTING] Server is starting.")
	client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	client.connect(ADDR)
	file = open("test.txt", "r")
	data = file.read()
	client.send("test.txt".encode(FORMAT))
	
	msg = client.recv(SIZE).decode(FORMAT)
	print("[Server] server replied")

	client.send(data.encode(FORMAT))
	msg = client.recv(SIZE).decode(FORMAT)
	print("[Server] server received data")
	
	file.close()
	client.close()

if __name__ == "__main__":
	main()
	

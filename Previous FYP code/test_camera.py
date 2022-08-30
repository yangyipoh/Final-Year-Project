import cv2

cap = cv2.VideoCapture(0)

if not cap.isOpened():
	raise IOError("Cannot open webcam")


while True:
	ret, frame = cap.read()
	if ret:
		cv2.imshow("Input", frame)
		if cv2.waitKey(1) == 27:
			break
	else:
		pass

cap.release()
cv2.destroyAllWindows()

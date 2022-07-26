import cv2
import numpy as np

cap = cv2.VideoCapture(1)
# cap.set(cv2.CAP_PROP_FRAME_WIDTH, 640)
# cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 480)
# cap.set(cv2.CAP_PROP_FPS, 30)

if (not cap.isOpened()):
	print("Unable to read camera feed")

fourcc = cv2.VideoWriter_fourcc("X", "V", "I", "D")
videoWriter = cv2.VideoWriter('video.avi', fourcc, 30, (640, 480))

while (True):
	ret, frame = cap.read()
	if ret:
		cv2.imshow('video', frame)
		videoWriter.write(frame)
		
	if cv2.waitKey(1) & 0xFF == ord('q'):
		break

cap.release()
videoWriter.release()

cv2.destroyAllWindows()

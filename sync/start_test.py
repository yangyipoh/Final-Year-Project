'''
IMPORTANT: Set your display resolution to 1920x1080
IMPORTANT: Open all the window in full screen

Reorder your task bar
Open CoolWin (baud rate @ 115200, connect it),
Open OpenBCI (write to SD, start session, configure electrodes, notch filter, BP filter)
'''

import pyautogui as auto
import cv2
import time

person_ID = 'A0'
session_ID = '0'

# grab the image
cap = cv2.VideoCapture(0)

# start recording
fourcc = cv2.VideoWriter_fourcc("X", "V", "I", "D")
videoWriter = cv2.VideoWriter('video.avi', fourcc, 30, (640, 480))

_ = input('Press Enter to start')

# Setup Arduino sampling
auto.hotkey('Win', '1')
auto.hotkey('ctrl', 'r')
time.sleep(0.2)
auto.write(f'{person_ID}_{session_ID}')
auto.press('Enter')

# setup OpenBCI
auto.hotkey('Win', '2')
auto.moveTo(92, 72)
auto.click()

# record
while (True):
	ret, frame = cap.read()
	if ret:
		cv2.imshow('video', frame)
		videoWriter.write(frame)
	if cv2.waitKey(1) & 0xFF == ord('q'):
		break

# stop video recording
cap.release()
videoWriter.release()
cv2.destroyAllWindows()

# # stop Arudino
auto.hotkey('Win', '3')
auto.hotkey('Win', '1')
auto.hotkey('ctrl', 'shift', 'r')

# stop OpenBCI
auto.hotkey('Win', '3')
auto.hotkey('Win', '2')
auto.moveTo(92, 72)
auto.click()





# ---------------------------------------------------------
# auto.moveTo(35, 100)
# auto.doubleClick()
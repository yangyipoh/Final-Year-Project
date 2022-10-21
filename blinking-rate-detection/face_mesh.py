"""
MediaPipe Code to extract facial features
From https://google.github.io/mediapipe/solutions/face_mesh
EAR https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9044337/
"""
import cv2
import mediapipe as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter,filtfilt, find_peaks
from tqdm import tqdm
import sys


LEFT_EYE = [362, 385, 387, 263, 373, 380]    # p1, p2, p3, p4, p5, p6
RIGHT_EYE = [33, 160, 158, 133, 153, 144]    # p1, p2, p3, p4, p5, p6
EAR_THRESHOLD = 0.020


def draw_features(img, landmark_lst, width, height, util):
    if not landmark_lst:
        return

    # plot features
    for idx in LEFT_EYE:
        landmark = landmark_lst.landmark[idx]
        img_coord = util._normalized_to_pixel_coordinates(landmark.x, landmark.y, width, height)
        cv2.circle(img, img_coord, 3, (255, 255, 255), -1)

    for idx in RIGHT_EYE:
        landmark = landmark_lst.landmark[idx]
        img_coord = util._normalized_to_pixel_coordinates(landmark.x, landmark.y, width, height)
        cv2.circle(img, img_coord, 3, (255, 255, 255), -1)
    
    return


def distance(point1, point2):
    x1, y1 = point1.x, point2.y
    x2, y2 = point2.x, point2.y
    return math.sqrt((x1-x2)**2 + (y1-y2)**2)


def count_blink(landmark_lst, verbose=False):
    if not landmark_lst:
        return 0
    
    # left eye
    vertical1_left = distance(landmark_lst.landmark[LEFT_EYE[1]], landmark_lst.landmark[LEFT_EYE[5]])
    vertical2_left = distance(landmark_lst.landmark[LEFT_EYE[2]], landmark_lst.landmark[LEFT_EYE[4]])
    horizontal_left = distance(landmark_lst.landmark[LEFT_EYE[0]], landmark_lst.landmark[LEFT_EYE[3]])

    # right eye
    vertical1_right = distance(landmark_lst.landmark[RIGHT_EYE[1]], landmark_lst.landmark[RIGHT_EYE[5]])
    vertical2_right = distance(landmark_lst.landmark[RIGHT_EYE[2]], landmark_lst.landmark[RIGHT_EYE[4]])
    horizontal_right = distance(landmark_lst.landmark[RIGHT_EYE[0]], landmark_lst.landmark[RIGHT_EYE[3]])

    ear_left = 0.5*(vertical1_left + vertical2_left)/horizontal_left
    # ear_right = 0.5*(vertical1_right + vertical2_right)/horizontal_right
    # ear_score = 0.5*(ear_left + ear_right)
    ear_score = ear_left

    # print(ear_score)
    if verbose and ear_score < EAR_THRESHOLD:
        print('Eye is closed')
    
    return ear_score

def butter_bandpass_filter(data, cutoff, fs, order):
    normal_cutoff = cutoff / (fs*0.5)
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='bandpass', analog=False)
    y = filtfilt(b, a, data)
    return y


def main():
    mp_drawing = mp.solutions.drawing_utils
    mp_face_mesh = mp.solutions.face_mesh
    # cap = cv2.VideoCapture(0)
    cap = cv2.VideoCapture(sys.argv[1])
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fs = cv2.CAP_PROP_FPS

    close_data = np.array([None]*total_frames)
    idx = 0
    with mp_face_mesh.FaceMesh(
        max_num_faces=1,
        refine_landmarks=True,
        min_detection_confidence=0.5,
        min_tracking_confidence=0.5) as face_mesh:
        with tqdm(total=total_frames) as pbar:
            while cap.isOpened():
                success, image = cap.read()
                if not success:
                    break

                # To improve performance, optionally mark the image as not writeable to
                # pass by reference.
                image.flags.writeable = False
                image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
                results = face_mesh.process(image)

                # Draw the face mesh annotations on the image.
                image.flags.writeable = True
                image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)
                if results.multi_face_landmarks:
                    face_landmarks = results.multi_face_landmarks[0]
                    # drawing eye on face
                    draw_features(image, face_landmarks, image.shape[1], image.shape[0], mp_drawing)

                    # count blinking rate
                    score = count_blink(face_landmarks, verbose=False)
                    close_data[idx] = score
                    idx += 1
                    pbar.update(1)
                    
                # cv2.imshow('MediaPipe Face Mesh', cv2.flip(image, 1))
                # if cv2.waitKey(5) & 0xFF == 27:
                #     break
                
            cap.release()
      
    # Create bandpass filter
    close_data = close_data[close_data != np.array(None)]
    filtered_data = butter_bandpass_filter(close_data, np.array([0.1, 0.7]), fs, 5)
    # filtered_data = close_data

    # find peaks
    peak_idx, _ = find_peaks(-filtered_data, height=0.01, distance=5, prominence=(None, None))
    print(f'Total blinks: {len(peak_idx)}')
    
    t_stamp = np.array([i/fs for i in range(len(filtered_data))])
    plt.plot(t_stamp, filtered_data, label='EAR')
    plt.plot(t_stamp[peak_idx], filtered_data[peak_idx], "x", label='Blink')
    plt.title(f'Total blinks: {len(peak_idx)}')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Eye Aspect Ratio')
    plt.legend()

    plt.show()

    # cv2.destroyAllWindows()


if __name__ == '__main__':
    main()
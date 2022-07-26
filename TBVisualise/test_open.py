from matplotlib import pyplot as plt
from scipy.ndimage import uniform_filter1d

f_wheel = open('Data/Wheel.txt')
f_acc = open('Data/Acc.txt')
f_brake = open('Data/Brake.txt')
f_time = open('Data/Time.txt')

# processing the data
N = 7
wheel_data = f_wheel.readlines()
wheel_data = [float(data.strip()) for data in wheel_data]
wheel_data = uniform_filter1d(wheel_data, size=N)
acc_data = f_acc.readlines()
acc_data = [float(data.strip()) for data in acc_data]
acc_data = uniform_filter1d(acc_data, size=N)
brake_data = f_brake.readlines()
brake_data = [float(data.strip()) for data in brake_data]
brake_data = uniform_filter1d(brake_data, size=N)
time_data = f_time.readlines()
time_data = [float(data.strip())/1e6 for data in time_data]

# plotting
# fig = plt.figure()
# ax1 = fig.add_subplot(311, constrained_layout=True)
# ax2 = fig.add_subplot(312, constrained_layout=True)
# ax3 = fig.add_subplot(313, constrained_layout=True)
# ax1.plot(time_data, acc_data)
# ax1.title.set_text('Accelerator Pedal data')
# ax2.plot(time_data, brake_data)
# ax2.title.set_text('Brake Pedal data')
# ax3.plot(time_data, wheel_data)
# ax3.title.set_text('Wheel data')
# plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
# plt.show()

fig, axs = plt.subplots(3, 1, constrained_layout=True)
ax1 = axs.flat[0]
ax2 = axs.flat[1]
ax3 = axs.flat[2]
ax1.plot(time_data, acc_data)
ax1.set_title('Accelerator Pedal data')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position')
ax2.plot(time_data, brake_data)
ax2.set_title('Brake Pedal data')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Position')
ax3.plot(time_data, wheel_data)
ax3.set_title('Wheel data')
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Steering angle')

# plt.plot(time_data, acc_data)
# plt.suptitle('Acceleration data')
# plt.subplot(3, 1, constrained_layout=True)
# plt.plot(time_data, brake_data)
# plt.subplot(3, 1, 3)
# plt.plot(time_data, wheel_data)
plt.show()

f_wheel.close()
f_acc.close()
f_brake.close()
f_time.close()

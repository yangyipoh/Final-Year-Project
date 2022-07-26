clear; close all; clc;

% Connect to serial port and set properties
portnum = serialportlist("available");
baud_rate = 115200;
serial1 = serialport(portnum, baud_rate);

size = 2000;
acc = zeros(size, 3);
idx = 1;
while idx <= size
    flush(serial1)
    recv_data = read(serial1, 20, 'string');
    acc(idx, :) = str2double(split(recv_data,',')) - 5;
    fprintf("Data: %s\n", recv_data)
    idx = idx + 1;
end

read_acc_data(acc);
clear variables; close all; clc;

blink_data = [367, 316, 270;
              487, 372, 206;
              496, 435, 240];

score = [5, 5.5, 4;
         7, 6, 5;
         8, 7, 7];

X = categorical({'Morning','Evening','Night'});
X = reordercats(X,{'Morning','Evening','Night'});
figure(1)
bar(X,blink_data)
xlabel('Time of day')
ylabel('Total Blinks')
title('Total blinks from three different times')
legend('Participant 1', 'Participant 2', 'Participant 3', 'location', 'nw')

blink_data = blink_data(:, 3);
score = score(:, 3);

y = blink_data(:);
x = score(:);

X = [ones(length(x),1) x];
b = X\y;

figure(2)
plot(score(:), blink_data(:), '*')
hold on
plot(x, X*b, '-')
xlabel('KSS score')
ylabel('Total blinks')
title('Participant 3')

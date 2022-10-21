import pygame
pygame.init()

# variables
WHITE = (255, 255, 255)
wheel_img = pygame.image.load('Asset/steering_wheel.jpg')
wheel_loc = (250, 250)
angle = 90

# open files
f_wheel = open('Data/Wheel.txt')
f_acc = open('Data/Acc.txt')
f_brake = open('Data/Brake.txt')
f_time = open('Data/Time.txt')

# Set up the drawing window
screen = pygame.display.set_mode([500, 500])

# Run until the user asks to quit
running = True
while running:

    # Did the user click the window close button?
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # Fill the background with white
    screen.fill(WHITE)

    # Draw circle
    rotated_image = pygame.transform.rotate(wheel_img, angle)
    new_rect = rotated_image.get_rect(center = wheel_img.get_rect(topleft = wheel_loc).center)

    screen.blit(rotated_image, new_rect)

    # Flip the display
    pygame.display.flip()

# Done! Time to quit.
pygame.quit()

# close files
f_wheel.close()
f_acc.close()
f_brake.close()
f_time.close()
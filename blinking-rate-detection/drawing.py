import sys
 
import pygame
from pygame.locals import *
 
pygame.init()
 
fps = 60
fpsClock = pygame.time.Clock()
 
width, height = 640, 480
screen = pygame.display.set_mode((width, height))
font = pygame.font.Font('freesansbold.ttf', 16)
text1 = font.render('P1', True, (0, 0, 0))
text2 = font.render('P2', True, (0, 0, 0))
text3 = font.render('P3', True, (0, 0, 0))
text4 = font.render('P4', True, (0, 0, 0))
text5 = font.render('P5', True, (0, 0, 0))
text6 = font.render('P6', True, (0, 0, 0))

# Game loop.
while True:
  screen.fill((255, 255, 255))
  
  for event in pygame.event.get():
    if event.type == QUIT:
      pygame.quit()
      sys.exit()
  
  # Update.
  
  # Draw.
  pygame.draw.ellipse(screen, (0, 0, 255), pygame.Rect(100, 200, 300, 50), 3)
  pygame.draw.ellipse(screen, (0, 0, 255), pygame.Rect(200, 200, 100, 50))
  pygame.draw.ellipse(screen, (255, 0, 0), pygame.Rect(180-12, 200-4, 24, 24)) #p2
  pygame.draw.ellipse(screen, (255, 0, 0), pygame.Rect(320-12, 200-4, 24, 24)) #p3
  pygame.draw.ellipse(screen, (255, 0, 0), pygame.Rect(180-12, 285-4-50, 24, 24)) #p6
  pygame.draw.ellipse(screen, (255, 0, 0), pygame.Rect(320-12, 285-4-50, 24, 24)) #p5
  pygame.draw.ellipse(screen, (255, 0, 0), pygame.Rect(100-12, 250-12-25, 24, 24)) #p1
  pygame.draw.ellipse(screen, (255, 0, 0), pygame.Rect(400-12, 250-12-25, 24, 24)) #p4
  pygame.draw.line(screen, (255, 0, 0), (180, 208), (180, 293-50), 6)
  pygame.draw.line(screen, (255, 0, 0), (320, 293-50), (320, 208), 6)
  pygame.draw.line(screen, (255, 0, 0), (100, 250-25), (400, 250-25), 6)
  screen.blit(text1, pygame.Rect(100-12+3, 250-12+3-25, 24, 24))
  screen.blit(text2, pygame.Rect(180-12+3, 200-4+3, 24+3, 24))
  screen.blit(text3, pygame.Rect(320-12+3, 200-4+3, 24, 24))
  screen.blit(text4, pygame.Rect(400-12+3, 250-12+3-25, 24, 24))
  screen.blit(text5, pygame.Rect(320-12+3, 285-4+3-50, 24, 24))
  screen.blit(text6, pygame.Rect(180-12+3, 285-4+3-50, 24, 24))
  
  pygame.display.flip()
  fpsClock.tick(fps)
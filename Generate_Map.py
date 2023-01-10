# !/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math
import cv2
import random


def Height_Matrix_Generate(points: int, size_matrix: int, grid_size: int, h: int, mode: int):
    # points - liczba punktów do tworzenia wzgórz
    # size_matrix - rozmiar macierzy -> mapy
    # grid_size - rozmiar siatki przy generowaniu macierzy
    # h - promień generowanych wzniesień
    # mode - tryb:
        # Mode 1 - tryb z wybieraniem punktów
        # Mode 0 - tryb losowego wyboru punktów

    # Function to calculate intensity with Quartic Kernel
    def kde_quartic(d, h):
        dn = d / h
        P = (15 / 16) * (1 - dn ** 2) ** 2
        return int(P * 10)

    # Create point matrix get coordinates of mouse click on image
    def mousePoints(event, x, y, flags, params):
        # global counter
        # Left button mouse click event opencv
        if event == cv2.EVENT_LBUTTONDOWN:
            point_matrix[counter] = x, y
            # counter = counter + 1

    scale = 5  # Przeskalowanie rozmiaru macierzy w celu łatwiejszego naniesienia punktów
    # Saving points
    point_matrix = np.zeros((points, 2), int)
    # Mode 1 - tryb z wybieraniem punktów
    # Mode 0 - tryb losowego wyboru punktów
    if mode == 1:
        scale = 5  # Przeskalowanie rozmiaru macierzy w celu łatwiejszego naniesienia punktów
        # Saving points
        point_matrix = np.zeros((points, 2), int)
        counter = 0

        # Read image
        img = np.zeros(shape=(size_matrix * scale, size_matrix * scale, 3))
        img = np.array(img, dtype=np.uint8)

        while counter < len(point_matrix):
            for x in range(0, points):
                cv2.circle(img, (point_matrix[x][0], point_matrix[x][1]), 3, (0, 255, 0), cv2.FILLED)
            # Showing original image
            cv2.imshow("Original Image ", img)
            # Mouse click event on original image
            cv2.setMouseCallback("Original Image ", mousePoints, point_matrix)
            # zliczenie dodanych punktów
            cnt = 0
            for i in point_matrix:
                if i[0] != 0 or i[1] != 0:
                    cnt += 1
            if cnt > counter:
                counter += 1
            # Refreshing window all time
            cv2.waitKey(1)
        cv2.destroyAllWindows()
        point_matrix = (point_matrix / scale).astype(int)
    elif mode == 0:
        for i in range(points):
            point_matrix[i, 0] = random.randint(0, size_matrix)
            point_matrix[i, 1] = random.randint(0, size_matrix)

    x = point_matrix[:, 0]
    y = point_matrix[:, 1]

    # Getting X,Y min and max
    x_min = 0
    x_max = size_matrix
    y_min = 0
    y_max = size_matrix

    # Construct grid
    x_grid = np.arange(x_min, x_max, grid_size)
    y_grid = np.arange(y_min, y_max, grid_size)
    x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)

    # Grid center point
    xc = x_mesh + (grid_size / 2)
    yc = y_mesh + (grid_size / 2)

    # Processing
    intensity_list = []
    for j in range(len(xc)):
        intensity_row = []
        for k in range(len(xc[0])):
            kde_value_list = []
            for i in range(len(x)):
                # Calculate distance
                d = math.sqrt((xc[j][k] - x[i]) ** 2 + (yc[j][k] - y[i]) ** 2)
                if d <= h:
                    p = kde_quartic(d, h)
                else:
                    p = 0
                kde_value_list.append(p)
            # sum all intensity value
            p_total = sum(kde_value_list)
            intensity_row.append(p_total)
        intensity_list.append(intensity_row)

    # HeatMap Output
    HeatMatrix = np.array(intensity_list)
    HeatMatrix = (200 * HeatMatrix / np.max(HeatMatrix)).astype(int)
    return HeatMatrix, x, y


def Terrain_Matrix_Generate(size_matrix, size_terrain_list, base_terrain, extra_terrains):
    # size_matrix - rozmiar macierzy -> mapy
    # size_terrain_list - możliwe rozmiary terenów
    # base_terrain - trudność terenu bazowego
    # extra_terrains - rodzaje terenu w postaci poziomów trudności
    terrainMatrix = np.full((size_matrix, size_matrix), base_terrain, dtype=int)
    for i in range(100):
        size = random.choice(size_terrain_list)  # Rozmiar wyspy danego terenu
        pos_x = random.randint(0, size_matrix - size + 1)
        pos_y = random.randint(0, size_matrix - size + 1)
        terrain = random.choice(extra_terrains)
        terrainMatrix[pos_x: pos_x + size, pos_y: pos_y + size] = terrain
    return terrainMatrix


def Generate_Matrices(size_matrix=100, points=30, grid_size=1, h=25, mode=0, size_terrain_list=[20, 25, 30, 35, 38],
                      base_terrain=10, extra_terrains=[10, 20, 30, 40]):
    HeatMatrix, x, y = Height_Matrix_Generate(points, size_matrix, grid_size, h, mode)
    terrainMatrix = Terrain_Matrix_Generate(size_matrix, size_terrain_list, base_terrain, extra_terrains)
    return HeatMatrix, terrainMatrix, x, y

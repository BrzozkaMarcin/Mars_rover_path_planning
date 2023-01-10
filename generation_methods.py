#!/usr/bin/python
from typing import List
import random
from collections import namedtuple

Point = namedtuple('Point', ('x', 'y'))


def Find_Random_Solution(FirstPoint, LastPoint):
    Solution = [FirstPoint]
    CurrentPoint = FirstPoint
    # Find pseudorandom solution
    while CurrentPoint.x != LastPoint.x or CurrentPoint.y != LastPoint.y:
        if CurrentPoint.x < LastPoint.x:
            if CurrentPoint.y < LastPoint.y:
                direction = random.randint(0, 2)
                if direction == 0:
                    CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y)
                elif direction == 1:
                    CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y + 1)
                elif direction == 2:
                    CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y + 1)
            elif CurrentPoint.y == LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y)
            elif CurrentPoint.y > LastPoint.y:
                direction = random.randint(0, 2)
                if direction == 0:
                    CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y)
                elif direction == 1:
                    CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y - 1)
                elif direction == 2:
                    CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y - 1)
        elif CurrentPoint.x == LastPoint.x:
            if CurrentPoint.y < LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y + 1)
            elif CurrentPoint.y > LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y - 1)
        elif CurrentPoint.x > LastPoint.x:
            if CurrentPoint.y < LastPoint.y:
                direction = random.randint(0, 2)
                if direction == 0:
                    CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y)
                elif direction == 1:
                    CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y + 1)
                elif direction == 2:
                    CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y + 1)
            elif CurrentPoint.y == LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y)
            elif CurrentPoint.y > LastPoint.y:
                direction = random.randint(0, 2)
                if direction == 0:
                    CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y)
                elif direction == 1:
                    CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y - 1)
                elif direction == 2:
                    CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y - 1)
        Solution.append(CurrentPoint)
    return Solution


def GenerateNeighbourSolution1(solution: List[Point], number_of_points: int = 10,
                               SizeMap: int = 10, bad_point_idx: int = 10, rng: int = 10) -> List[Point]:
    neighbour_solution = [solution[0]]
    # Select n random points on a path
    points = random.sample(solution, min(len(solution) - 1,number_of_points))
    static_points = []
    for el in solution:
        if el in points:
            static_points.append(el)
    static_points = [solution[0]] + static_points + [solution[-1]]
    # For every pair of points, generate a random path between them
    for i in range(len(static_points) - 1):
        partial_solution = Find_Random_Solution(static_points[i], static_points[i + 1])
        neighbour_solution += partial_solution[1:]
    return neighbour_solution


def GenerateNeighbourSolution2(solution: List[Point], number_of_points: int = 10,
                               SizeMap: int = 10, bad_point_idx: int = 10, rng: int = 10) -> List[Point]:
    half_rng = int(rng / 2)
    id_start = max(bad_point_idx - half_rng, 1)
    id_end = min(bad_point_idx + half_rng, len(solution) - 2)
    StartPoint = solution[id_start]
    EndPoint = solution[id_end]
    S1 = solution[:id_start + 1]
    S4 = solution[id_end:]
    idx, idy = solution[bad_point_idx].x, solution[bad_point_idx].y
    idx_min, idx_max = max(idx - half_rng, 0), min(idx + half_rng, SizeMap - 1)
    if idx_min > idx_max: idx_min, idx_max = idx_max, idx_min
    idy_min, idy_max = max(idy - half_rng, 0), min(idy + half_rng, SizeMap - 1)
    if idy_min > idy_max: idy_min, idy_max = idy_max, idy_min
    IntermediatePoint = Point(random.randint(
        idx_min, idx_max), random.randint(idy_min, idy_max))
    static_points = [solution[0], StartPoint,
                     IntermediatePoint, EndPoint, solution[-1]]
    S2 = Find_Random_Solution(StartPoint, IntermediatePoint)
    S3 = Find_Random_Solution(IntermediatePoint, EndPoint)
    neighbour_solution = S1 + S2[1:] + S3[1:] + S4[1:]
    # return neighbour_solution, static_points
    return neighbour_solution


def GenerateNeighbourSolution3(solution: List[Point], number_of_points: int = 10,
                               SizeMap: int = 10, bad_point_idx: int = 10, rng: int = 10) -> List[Point]:
    half_rng = int(rng / 2)
    id_start = max(bad_point_idx - half_rng, 1)
    id_end = min(bad_point_idx + half_rng, len(solution)-2)
    StartPoint = solution[id_start]
    LastPoint = solution[id_end]
    S1 = solution[:id_start + 1]
    S3 = solution[id_end:]
    static_points = [solution[0], StartPoint, LastPoint, solution[-1]]
    S2 = [StartPoint]
    CurrentPoint = StartPoint
    state = random.choice([False, True])
    if not state:
        while CurrentPoint.x != LastPoint.x:
            if CurrentPoint.x < LastPoint.x:
                CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y)
            elif CurrentPoint.x > LastPoint.x:
                CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y)
            S2.append(CurrentPoint)
        while CurrentPoint.y != LastPoint.y:
            if CurrentPoint.y < LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y + 1)
            elif CurrentPoint.y > LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y - 1)
            S2.append(CurrentPoint)
    else:
        while CurrentPoint.y != LastPoint.y:
            if CurrentPoint.y < LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y + 1)
            elif CurrentPoint.y > LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y - 1)
            S2.append(CurrentPoint)
        while CurrentPoint.x != LastPoint.x:
            if CurrentPoint.x < LastPoint.x:
                CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y)
            elif CurrentPoint.x > LastPoint.x:
                CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y)
            S2.append(CurrentPoint)
    neighbour_solution = S1 + S2[1:] + S3[1:]
        # return neighbour_solution, static_points
    return neighbour_solution


def GenerateNeighbourSolution4(solution: List[Point], number_of_points: int = 10,
                               SizeMap: int = 10, bad_point_idx: int = 10, rng: int = 20) -> List[Point]:
    half_rng = int(rng / 2)
    id_start = max(bad_point_idx - half_rng, 1)
    id_end = min(bad_point_idx + half_rng, len(solution)-2)
    StartPoint = solution[id_start]
    LastPoint = solution[id_end]
    S1 = solution[:id_start + 1]
    S3 = solution[id_end:]
    static_points = [solution[0], StartPoint, LastPoint, solution[-1]]
    S2 = [StartPoint]
    CurrentPoint = StartPoint
    state = random.choice([False, True])
    while CurrentPoint.x != LastPoint.x or CurrentPoint.y != LastPoint.y:
        if not state:
            if CurrentPoint.x < LastPoint.x:
                CurrentPoint = Point(CurrentPoint.x + 1, CurrentPoint.y)
            elif CurrentPoint.x > LastPoint.x:
                CurrentPoint = Point(CurrentPoint.x - 1, CurrentPoint.y)
            state = True
        else:
            if CurrentPoint.y < LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y + 1)
            elif CurrentPoint.y > LastPoint.y:
                CurrentPoint = Point(CurrentPoint.x, CurrentPoint.y - 1)
            state = False
        S2.append(CurrentPoint)
    neighbour_solution = S1 + S2[1:] + S3[1:]
    # if neighbour_solution == solution:
    #     neighbour_solution, static_points = GenerateNeighbourSolution4(
    #         solution, bad_point_idx, rng)
    # return neighbour_solution, static_points
    return neighbour_solution

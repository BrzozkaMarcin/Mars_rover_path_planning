#!/usr/bin/python
# -*- coding: utf-8 -*-
from math import inf
from generation_methods import *

Point = namedtuple('Point', ('x', 'y'))


# Class for solution
class SolutionClass:
    def __init__(self, Solution=None, time=None):
        self.Solution = Solution
        self.time = time


# Class with algorithm
class TabuSearch:
    def __init__(self, HeightMatrix, TerrainMatrix, StartPoint, EndPoint, MaxEnergy, iteration,
                 height_constant_energy, terrain_constant_energy, height_constant_time, terrain_constant_time,
                 NeighbourMethodList=None, sn1=8, sn2=4, aspiration=0, randompoint=5, neighbourhood_size=20,
                 k=30, k_p=10, war_zakonczenia=0):
        if NeighbourMethodList is None:
            self.GenerationMethodsNeighbourList = [GenerateNeighbourSolution1, GenerateNeighbourSolution2,
                                                   GenerateNeighbourSolution3, GenerateNeighbourSolution4]
        self.HeightMatrix = HeightMatrix
        self.SizeMap = len(HeightMatrix)
        self.TerrainMatrix = TerrainMatrix
        self.StartPoint = StartPoint
        self.EndPoint = EndPoint
        self.MaxEnergy = MaxEnergy
        self.height_constant_energy = height_constant_energy
        self.terrain_constant_energy = terrain_constant_energy
        self.height_constant_time = height_constant_time
        self.terrain_constant_time = terrain_constant_time
        self.FirstSolution = None
        self.BestSolution = None
        self.ActualSolution = None
        self.NeighbourList = []
        self.TabuList = []
        self.TabuListPoionts = []
        self.iteration = iteration
        self.GenerationMethodsNeighbourList = NeighbourMethodList
        self.GoalFunctionValues = []
        self.BestSolutionValues = []
        self.sn1 = sn1
        self.sn2 = sn2
        if aspiration == 0:
            self.aspiration = inf
        else:
            self.aspiration = aspiration
        self.random_point = randompoint
        self.neighbourhood_size = neighbourhood_size
        self.k = k
        self.k_p = k_p
        if war_zakonczenia == 0:
            self.war_zakonczenia = inf
        else:
            self.war_zakonczenia = war_zakonczenia

    # Calculate the time it takes to travel through one point
    def calculate_time(self, LastPoint, NewPoint):
        height_last = self.HeightMatrix[LastPoint.x, LastPoint.y]
        height_new = self.HeightMatrix[NewPoint.x, NewPoint.y]

        terrain_last = self.TerrainMatrix[LastPoint.x, LastPoint.y]
        terrain_new = self.TerrainMatrix[NewPoint.x, NewPoint.y]

        # If the cart drives downhill, change this constant
        if height_last > height_new:
            height_constant = 0.6 * self.height_constant_time  # Speed x2 if it goes down and time 2x shorter
        else:
            height_constant = self.height_constant_time

        height_time = height_constant * abs(height_new - height_last)
        terrain_time = self.terrain_constant_time * (terrain_new + terrain_last) / 2

        return height_time + terrain_time

    def calculate_total_time(self, Solution):
        TotalTime = 0
        for idx in range(1, len(Solution)):
            LastPoint = Solution[idx - 1]
            NewPoint = Solution[idx]
            TotalTime += self.calculate_time(LastPoint, NewPoint)
        return TotalTime

    # Calculate the energy change when driving through a point
    def calculate_energy(self, LastPoint, NewPoint):
        height_last = self.HeightMatrix[LastPoint.x, LastPoint.y]
        height_new = self.HeightMatrix[NewPoint.x, NewPoint.y]

        terrain_last = self.TerrainMatrix[LastPoint.x, LastPoint.y]
        terrain_new = self.TerrainMatrix[NewPoint.x, NewPoint.y]

        # If the cart drives downhill, change this constant
        if height_last > height_new:
            height_constant = -0.5 * self.height_constant_energy  # energy charging by 0.5
        else:
            height_constant = self.height_constant_energy

        height_energy = height_constant * abs(height_new - height_last)
        terrain_energy = self.terrain_constant_energy * (terrain_new + terrain_last) / 2

        # energy from sun (if x and y new point is > : 2, if x or y new point > : 1, else 0
        if LastPoint.x < NewPoint.x and LastPoint.y < NewPoint.y:
            energy_sun = self.sn1
        elif LastPoint.x < NewPoint.x or LastPoint.y < NewPoint.y:
            energy_sun = self.sn2
        else:

            energy_sun = 0
        return terrain_energy + height_energy - energy_sun

    # A function that checks whether the rover for this route will discharge before reaching the end point.
    def try_energy(self, Solution):
        energy = self.MaxEnergy
        cnt = 0
        for idx in range(1, len(Solution)):
            LastPoint = Solution[idx - 1]
            NewPoint = Solution[idx]
            calc_energy = self.calculate_energy(LastPoint, NewPoint)
            energy -= calc_energy
            cnt += calc_energy
            if energy <= 0:
                return 1  # Solution does not meet energy requirements
        return 0  # Solution is ok

    def first_solution(self):
        max_it = 100  # max number of try find first solution
        it = 0
        while it < max_it:
            Solution = [self.StartPoint]
            points = [Point(random.randint(1, 99), random.randint(1, 99)) for _ in range(2)]
            points = sorted(points, key=lambda k: [k[0], k[1]])
            points = [self.StartPoint] + points + [self.EndPoint]
            for i in range(len(points) - 1):
                partial_solution = Find_Random_Solution(points[i], points[i + 1])
                Solution += partial_solution[1:]
            TotalTime = self.calculate_total_time(Solution)
            if self.try_energy(Solution):
                it += 1
                continue
            else:
                return SolutionClass(Solution, TotalTime)
        return None

    # Algorithm TabuSearch
    def Tabu_Search_Algorithm(self):
        self.ActualSolution = self.first_solution()
        self.FirstSolution = self.ActualSolution
        if self.ActualSolution is None:
            print("Zwiększ energię pojazdu")
            return None, None

        counter_to_random_point = 0
        counter_to_aspiration = 0
        self.BestSolution = self.ActualSolution
        self.TabuList = []
        iterations_without_improvement = 0
        for _ in range(self.iteration):
            # Every a certain number of iterations look for a neighbor around a random point
            if counter_to_random_point == self.random_point:
                worst_point = random.randint(1, len(self.ActualSolution.Solution) - 2)
                counter_to_random_point = 0
            else:
                worst_point = self.WorstPoint(self.ActualSolution.Solution)
                counter_to_random_point += 1
            # Size hardcoded to 20 ,might be neccessary to change
            self.NeighbourList = self.GenerateNeighbourhood(base_solution=self.ActualSolution.Solution,
                                                            neighbourhood_size=self.neighbourhood_size,
                                                            bad_point_idx=worst_point)
            iterations_without_improvement += 1
            if iterations_without_improvement > self.war_zakonczenia:
                break
            BestNeighbour = SolutionClass(None, inf)
            for el in self.NeighbourList:
                if el in [x[0] for x in self.TabuList]:
                    continue
                if el.time < BestNeighbour.time:
                    BestNeighbour = el
            if BestNeighbour.Solution is None and self.TabuList:
                el_to_del = self.TabuList[0]
                for el in self.TabuList:
                    if el[0].time < BestNeighbour.time:
                        BestNeighbour = el[0]
                        el_to_del = el
                self.TabuList.remove(el_to_del)
                counter_to_aspiration = 0
            elif BestNeighbour.Solution is None and not self.TabuList:
                BestNeighbour = self.ActualSolution

            if BestNeighbour.time < self.BestSolution.time:
                self.BestSolution = BestNeighbour
                iterations_without_improvement = 0

            else:
                counter_to_aspiration += 1
                if counter_to_aspiration == self.aspiration:
                    BestNeighbour = SolutionClass(None, inf)
                    el_to_del = self.TabuList[0]
                    for el in self.TabuList:
                        if el[0].time < BestNeighbour.time:
                            BestNeighbour = el[0]
                            el_to_del = el
                    self.TabuList.remove(el_to_del)
                    counter_to_aspiration = 0

            self.ActualSolution = BestNeighbour
            # Dodawanie losowych punktów rozwiązania do listy punktów tabu
            # print(BestNeighbour.Solution)
            # print(max(len(BestNeighbour.Solution), 10))
            sample_points = random.sample(BestNeighbour.Solution, min(len(BestNeighbour.Solution) - 1, 10))
            k_p = self.k_p
            self.TabuListPoionts = self.TabuListPoionts + [[point, k_p] for point in sample_points]
            k = self.k
            self.TabuList.append([BestNeighbour, k])
            for Tabu in self.TabuList:
                Tabu[1] -= 1
                if Tabu[1] == 0:
                    self.TabuList.remove(Tabu)

            for Tabu in self.TabuListPoionts:
                Tabu[1] -= 1
                if Tabu[1] == 0:
                    self.TabuListPoionts.remove(Tabu)
            self.GoalFunctionValues.append(self.ActualSolution.time)
            self.BestSolutionValues.append(self.BestSolution.time)
        return self.FirstSolution, self.BestSolution

    # Generate a neighbourhood for a solution, of size n using selected method
    def GenerateNeighbourhood(self, base_solution: List[Point], neighbourhood_size: int,
                              number_of_points: int = 10, bad_point_idx: int = 30, rng: int = 10) -> List[
        SolutionClass]:
        neighbourhood = list()
        for _ in range(neighbourhood_size):
            generation_method = random.choice(self.GenerationMethodsNeighbourList)
            new_solution = generation_method(base_solution, number_of_points, self.SizeMap, bad_point_idx, rng)
            num_of_iterations = 50  # max number of try find neighbour
            add_solution = True
            while self.try_energy(new_solution) and not any(point in new_solution for point, _ in self.TabuListPoionts):
                new_solution = generation_method(base_solution)
                num_of_iterations -= 1
                if num_of_iterations <= 0:
                    add_solution = False
                    break
            if add_solution:
                new_solution_class = SolutionClass(new_solution, self.calculate_total_time(new_solution))
                neighbourhood.append(new_solution_class)
        if neighbourhood:
            return neighbourhood
        else:
            return []

    def WorstPoint(self, Solution) -> int:
        worst_Point = 0
        idx_worst_Point = None
        for idx in range(2, len(Solution) - 1):
            LastPoint = Solution[idx - 1]
            NewPoint = Solution[idx]
            point = self.calculate_time(LastPoint, NewPoint)
            if point > worst_Point:
                worst_Point = point
                idx_worst_Point = idx
        return idx_worst_Point

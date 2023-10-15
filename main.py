import pandas as pd
import os
import time
from tqdm import tqdm
import math
import matplotlib.pyplot as plt
import numpy as np
from pandas import *

FILES = ["09"]

SREZ = -1
data_coordinates = {}
data = {}

for i in FILES:
    print("Reading M" + i)
    file_path_coordinates = "data/HACK_OCT_2023/Coordinates/PositionsM" + i + ".dat"

    # Чтение файла и создание датафрейма
    df = pd.read_csv(file_path_coordinates, delimiter="\t")
    data_coordinates["M" + i] = df[0:SREZ:1]

    data["M" + i] = {}

    for file_name in tqdm(os.listdir("data/HACK_OCT_2023/M" + i)[0:SREZ:1], desc="Обработка файлов", ncols=100, ascii=True):
        df_event = pd.read_csv("data/HACK_OCT_2023/M" + i + "/" + file_name, delimiter="\t")
        data["M" + i][file_name.replace(".dat", "")] = df_event


class Particle:
    def __init__(self, energy_, x, y):
        self.energy = energy_
        self.xy = (x, y)


    def get_enery(self):
        return self.energy


    def get_position(self):
        return self.xy


class Detector:
    def __init__(self, size_of_unit_):
        self.size_of_unit = size_of_unit_
        self.height = 1000.0
        self.width = 1000.0
        self.count_of_particles = 0
        self.matrix = []
        for ii in range(int(self.height // self.size_of_unit)):
            mas = []
            for jj in range(int(self.width // self.size_of_unit)):
                mas.append([])
            self.matrix.append(mas.copy())

        self.matrix_part_counts = []


        for ii in range(int(self.height // self.size_of_unit)):
            mas = []
            for jj in range(int(self.width // self.size_of_unit)):
                mas.append(0)
            self.matrix_part_counts.append(mas)


    def get_matrix_energy(self):
        matrix_part_energies = []
        for ii in range(int(self.height // self.size_of_unit)):
            mas = []
            for jj in range(int(self.width // self.size_of_unit)):
                mas.append(self.matrix_part_counts[ii][jj] * self.get_mean_energy_of_particles_in_unit(jj, ii))
                # mas.append(self.get_mean_energy_of_particles_in_unit(ii, jj))
            matrix_part_energies.append(mas)
        return matrix_part_energies

    def add_particle(self, particle):
        new_x, new_y = particle.get_position()[0] + 500, particle.get_position()[1] + 500
        row, line = int(new_x // self.size_of_unit), int(new_y // self.size_of_unit)
        try:
            # Прямое добавление
            self.matrix[line][row].append(particle)
            self.matrix_part_counts[line][row] += 1
            self.count_of_particles += 1

            # Обратное добавление
            if (len(self.matrix) % 2 != 0):
                if (row == (len(self.matrix) // 2)):
                    pass
                else:
                    self.matrix[line][len(self.matrix) - row - 1].append(particle)
                    self.matrix_part_counts[line][len(self.matrix) - row - 1] += 1
                    self.count_of_particles += 1
            else:
                self.matrix[line][len(self.matrix) - row - 1].append(particle)
                self.matrix_part_counts[line][len(self.matrix) - row - 1] += 1
                self.count_of_particles += 1


        except:
            pass


    def get_static_weight_of_unit(self, i, j):
        return self.matrix_part_counts[j][i] / self.count_of_particles

    def get_mean_energy_of_particles_in_unit(self, i, j):
        sum_of_energy = 0
        for p in self.matrix[j][i]:
            sum_of_energy += p.get_enery()

        try:
            return sum_of_energy / self.matrix_part_counts[j][i]
        except:
            return 0

    def get_mean_energy_of_particles(self):
        sum_of_energy = 0
        for s in self.matrix:
            for stolb in s:
                for el in stolb:
                    sum_of_energy += el.get_enery()

        return sum_of_energy / self.count_of_particles

    def get_non_homogeneity(self):
        Q = self.get_mean_energy_of_particles()

        a = 0
        for st in range(0, len(self.matrix)):
            for stolb in range(0, len(self.matrix[0])):
                Wij = self.get_static_weight_of_unit(stolb, st)
                Qij = self.get_mean_energy_of_particles_in_unit(stolb, st)
                a += Wij * ((Qij - Q) ** 2)


        non_homogeneity = (math.sqrt(a) / Q) * 100

        # Нарисовать матрицу светосбора
        plt.imshow(self.get_matrix_energy(), cmap='viridis')  # cmap - цветовая карта (может быть изменена)
        plt.colorbar()  # Добавление цветовой шкалы

        plt.show()

        return non_homogeneity


detector = Detector(20)
for file in FILES:
    for coordinates in tqdm(data_coordinates["M" + file].itertuples(), desc="Обработка файлов", ncols=100, ascii=True):
        x = coordinates[1]
        y = coordinates[2]
        try:
            dd = str(int(file + "00000") + int(coordinates[0]))
            if len(dd) < 7:
                dd = "0" + dd
            dd = "EventM" + dd
            if (file == "24" or file == "25" or file == "26"):
                dd = "M" + file + "Event" + "0" * (5 - len(str(coordinates[0]))) + str(coordinates[0])

            for e in data["M" + file][dd]["Energy"][:SREZ:]:
                particle = Particle(e, x, y)
                detector.add_particle(particle)
        except:
            print("ERROR", dd)

print("Коэффициент неоднородности: ", detector.get_non_homogeneity())

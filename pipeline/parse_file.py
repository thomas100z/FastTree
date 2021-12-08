from __future__ import annotations
import numpy as np


class Profile:

    def __init__(self, alignment: Alignment) -> None:
        self.profile = Profile.form_profile([alignment.alignment])

    @staticmethod
    def form_profile(sequences: list[str], psuesdocount=False) -> np.array[float]:
        base_value = 4 if psuesdocount else 0

        bases = ['A', 'C', 'G', 'T']

        # start with one due to pseudocounts
        if psuesdocount:
            profile = np.ones(shape=(4, len(sequences[0])))
        else:
            profile = np.zeros(shape=(4, len(sequences[0])))

        # calculate counts
        for i in range(0, len(sequences)):
            for j in range(0, len(sequences[i])):

                # if found gap, continue, can change so the probabilities sum up to 1
                if sequences[i][j] == " " or sequences[i][j] == "-":
                    continue
                profile[bases.index(sequences[i][j])][j] += 1

        return profile / (len(sequences) + base_value)

    def join_profiles(self, other_profile: Profile) -> Profile:
        return np.mean(np.array([self.profile, other_profile.profile]), axis=0)

    def distance(self, other_profile):
        S = 0
        for i in range(4):
            for j in range(len(self.profile[0])):
                S += abs(self.profile[i, j] - other_profile[i, j])

        return S / 4


class Alignment:

    def __init__(self, name: str, alignment: str) -> None:
        self.alignment = alignment
        self.name = name
        self.profile = Profile(self)


class AlignmentParser:

    def __init__(self, file_name: str) -> None:
        self.sequences = []
        lines = open(file_name).read().splitlines()

        for i in range(round(len(lines) / 2)):
            self.sequences.append(Alignment(lines[(i * 2)].strip(), lines[(i * 2) + 1].strip()))

    def get_data(self) -> list[Alignment]:
        return self.sequences

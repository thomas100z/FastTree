import math


class DistanceProfile:

    @staticmethod
    def uncorrected_distance(sequence_1: str, sequence_2: str) -> float:
        result = 0
        gaps = 0
        for x, y in zip(sequence_1, sequence_2):
            if x != '-' and y != '-':
                result += (x != y)
            else:
                gaps += 1

        return result / (len(sequence_1) - gaps)

    @staticmethod
    def corrected_distance(sequence_1: str, sequence_2: str) -> float:
        du = 1 - (4 / 3) * DistanceProfile.uncorrected_distance(sequence_1, sequence_2)
        if 0 < du:
            return -3 / 4 * math.log10(du)
        else:
            # TODO:
            # ask
            return 1.5

    @staticmethod
    def out_distance(sequence_1: str, sequence_2: str) -> float:
        return sum([x != y for x, y in zip(sequence_1, sequence_2)]) / (len(sequence_1) - 2)

    @staticmethod
    def up_distance():
        pass

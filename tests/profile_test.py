import unittest
from pipeline import Profile, Alignment, DistanceProfile


class ProfileTest(unittest.TestCase):

    def test_corrected_distance(self):
        A = Alignment('A', 'CA')
        B = Alignment('B', 'GG')
        C = Alignment('C', 'CC')

        ac = C.profile.distance(A.profile.profile)
        ab = B.profile.distance(A.profile.profile)
        a = C.profile.distance(A.profile.join_profiles(B.profile))
        ac_seq_dis = DistanceProfile.uncorrected_distance(A.alignment, C.alignment)
        ab_seq_dis = DistanceProfile.uncorrected_distance(A.alignment, B.alignment)


if __name__ == '__main__':
    unittest.main()

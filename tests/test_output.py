import sys

sys.path.insert(0, "..")
import os
import unittest
from vic.vic import volume
from vic.vic.utils import niftireader

OUTPUT_FILE = 'data/_output/Grouped.csv'


def _get_file_name(case_id, lesion_id, type):
    return os.path.join('tests', 'test_data', case_id, lesion_id,
                        '{0}_L{1}_{2}.nii.gz'.format(case_id, lesion_id, type))


class TestMargins(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)

    def test_01_no_overlap_10mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '01_no_overlap_10mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '01_no_overlap_10mm_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]
        self.assertEqual(record["VIC_0"], 0.515)
        self.assertEqual(record["VIC_5"], 3.809)
        self.assertEqual(record["VIC_10"], 13.546)


    def test_02_perfect_overlap_0mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '02_perfect_overlap_0mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '02_perfect_overlap_0mm_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 3.29, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 13.03, delta=0.01)

    def test_03_perfect_overlap_10mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '03_perfect_overlap_10mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '03_perfect_overlap_10mm_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 0, delta=0.01)

    def test_04_perfect_overlap_5mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '04_perfect_overlap_5mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '04_perfect_overlap_5mm_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 9.378, delta=0.01)

    def test_05_perfect_overlap_n5mm_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '05_perfect_overlap_-5mm_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '05_perfect_overlap_-5mm_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 3.654, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 13.032, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 31.878, delta=0.01)

    def test_06_perfect_shifted_5mm_xy_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '06_perfect_shifted_5mm_xy_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '06_perfect_shifted_5mm_xy_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.087, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 1.881, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 9.698, delta=0.01)

    def test_07_perfect_shifted_5mm_x_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '07_perfect_shifted_5mm_x_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '07_perfect_shifted_5mm_x_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 1.309, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 9.378, delta=0.01)

    def test_08_perfect_shifted_5mm_y_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '08_perfect_shifted_5mm_y_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '08_perfect_shifted_5mm_y_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 1.309, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 9.378, delta=0.01)

    def test_09_perfect_shifted_5mm_z_margin(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T01', '09_perfect_shifted_5mm_z_margin', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T01', '09_perfect_shifted_5mm_z_margin', 'Ablation'))
        vic_ml = volume.compute_vic(tumor_np, ablation_np, None, [1, 1, 1])
        df = volume.summarize_vic('T01', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T01") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 1.309, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 9.378, delta=0.01)


class TestSubcapsularLesions(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)

    def test_01_0mm_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '01_0mm_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '01_0mm_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '01_0mm_margin_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 1.489, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 5.574, delta=0.01)

    def test_02_5mm_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '02_5mm_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '02_5mm_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '02_5mm_margin_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0.468, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 4.537, delta=0.01)

    def test_03_n5mm_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '03_-5mm_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '03_-5mm_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '03_-5mm_margin_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 1.037, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 4.446, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 11.071, delta=0.01)

    def test_04_2mm_shifted_tumor_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '04_2mm_shifted_tumor_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '04_2mm_shifted_tumor_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '04_2mm_shifted_tumor_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0.764, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 5.569, delta=0.01)

    def test_05_4mm_shifted_tumor_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '05_4mm_shifted_tumor_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '05_4mm_shifted_tumor_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '05_4mm_shifted_tumor_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.0, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 1.132, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 6.757, delta=0.01)

    def test_06_shifted_ablation_5mm_xy_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '06_shifted_ablation_5mm_xy_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '06_shifted_ablation_5mm_xy_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '06_shifted_ablation_5mm_xy_margin_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.061, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0.834, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 4.214, delta=0.01)

    def test_07_shifted_5mm_x_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '07_shifted_5mm_x_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '07_shifted_5mm_x_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '07_shifted_5mm_x_margin_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.024, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0.708, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 3.992, delta=0.01)

    def test_08_shifted_5mm_y_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '08_shifted_ablation_5mm_y_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '08_shifted_ablation_5mm_y_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '08_shifted_ablation_5mm_y_margin_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.041, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0.867, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 4.587, delta=0.01)

    def test_09_shifted_ablation_5mm_z_margin_subcapsular(self):
        _, tumor_np = niftireader.load_image(_get_file_name('T02', '09_shifted_ablation_5mm_z_margin_subcapsular', 'Tumor'))
        _, ablation_np = niftireader.load_image(_get_file_name('T02', '09_shifted_ablation_5mm_z_margin_subcapsular', 'Ablation'))
        _, liver_np = niftireader.load_image(_get_file_name('T02', '09_shifted_ablation_5mm_z_margin_subcapsular', 'Liver'))

        vic_ml = volume.compute_vic(tumor_np, ablation_np, liver_np, [1, 1, 1])
        df = volume.summarize_vic('T02', 'L01', vic_ml)
        print(df)

        record = df.loc[(df['Subject'] == "T02") & (df['Tumor'] == 'L01')].iloc[0]

        self.assertAlmostEqual(record["VIC_0"], 0.041, delta=0.01)
        self.assertAlmostEqual(record["VIC_5"], 0.867, delta=0.01)
        self.assertAlmostEqual(record["VIC_10"], 4.587, delta=0.01)


if __name__ == '__main__':
    unittest.main()

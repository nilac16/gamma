import pydicom as pyd
import numpy as np
import matplotlib.pyplot as plt
import gamma
import math
import time


def load_distribution(ds: pyd.Dataset) -> gamma.Distribution:
    res = gamma.Distribution()

    e0 = np.array(ds.ImageOrientationPatient[0:3])
    e1 = np.array(ds.ImageOrientationPatient[3:6])
    e2 = np.cross(e0, e1)
    res.matrix = np.column_stack((e0, e1, e2))

    res.origin = np.array([x for x in ds.ImagePositionPatient])

    res.spacing = np.array([
        ds.PixelSpacing[0],
        ds.PixelSpacing[1],
        ds.SliceThickness
    ])

    scal = float(ds.DoseGridScaling)
    res.data = np.array(ds.pixel_array, dtype=np.double) * scal
    res.data = res.data.reshape((
        ds.Columns,
        ds.Rows,
        ds.NumberOfFrames
    ))
    return res


if __name__ == "__main__":
    import sys

    params = gamma.Parameters(
        diff      = 0.03,
        dta       = 2.0,
        threshold = 0.10,
        relative  = False
    )
    opts = gamma.Options(pattern_shrinks = 6)
    try:
        tps = load_distribution(pyd.dcmread(sys.argv[1]))
        meas = load_distribution(pyd.dcmread(sys.argv[2]))
    except IndexError:
        print("Supply a TPS and measurement DICOM")
        exit(1)

    #print(f"TPS dose max:      {np.max(tps.data)}")
    #print(f"{tps.matrix=}\n{tps.origin=}\n{tps.spacing=}\n{tps.data.shape}")
    #print(f"Measured dose max: {np.max(meas.data)}")
    #print(f"{meas.matrix=}\n{meas.origin=}\n{meas.spacing=}\n{meas.data.shape}")

    t0 = time.time()
    res = gamma.compute(params, opts, tps, meas)
    t1 = time.time()
    print(f"{'RELATIVE DOSE' if params.relative else 'GLOBAL DOSE'}\n"
          f"Passed:         {res.passed}\n"
          f"Total:          {res.total}\n"
          f"Rate:           {100.0 * float(res.passed) / float(res.total):.2f}%\n"
          f"Minimum:        {res.min}\n"
          f"Maximum:        {res.max}\n"
          f"Average:        {res.mean}\n"
          f"Std. deviation: {math.sqrt(res.msqr) - res.mean * res.mean}")
    print(f"Took {t1 - t0:.2f} s")

    plt.pcolor(res.dist[:, 50, :], shading="flat")
    plt.colorbar()
    plt.show()

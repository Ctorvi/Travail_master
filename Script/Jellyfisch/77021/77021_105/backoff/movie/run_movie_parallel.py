import numpy as np
from pathlib import Path
import sys
import multiprocessing
import subprocess





def call_make_movie_frame(angle, frame_index):
    print(f"Making movie frame for angle {angle}")
    subprocess.run(["python", "./Script/Jellyfisch/77021/77021_105/backoff/movie/make_movie_frame_JF.py", str(angle), str(frame_index)])


N_PROCS = 6
MAX_ANGLE = np.pi
NUM_FRAMES = 8


# if __name__ == "__main__":
#     angles = np.linspace(0., MAX_ANGLE, NUM_FRAMES)
#     with multiprocessing.Pool(4) as pool:
#         pool.starmap(call_make_movie_frame, [(angle, frame_index) for frame_index, angle in enumerate(angles)])
#     print("All frames done")


if __name__ == "__main__":
    angles = np.linspace(0., MAX_ANGLE, NUM_FRAMES)
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        pool.starmap(call_make_movie_frame, [(angle, frame_index) for frame_index, angle in enumerate(angles)])
    print("All frames done")
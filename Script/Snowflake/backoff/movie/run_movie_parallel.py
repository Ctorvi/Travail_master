import numpy as np
from pathlib import Path
import sys
import multiprocessing
import subprocess





def call_make_movie_frame(angle, frame_index):
    print(f"Making movie frame for angle {angle}")
    subprocess.run(["python", "./Script/Snowflake/movie/make_movie_frame_2.py", str(angle), str(frame_index)])


N_PROCS = 6
MAX_ANGLE = 2*np.pi
NUM_FRAMES = 1


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
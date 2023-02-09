from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt

import mhm

here = Path(__file__).parent

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
first_time = True

# disable terminal output of mHM
mhm.model.set_verbosity(level=1)
# use single test domain 1
mhm.model.init(cwd=here / ".." / ".." / "test_domain")
# disable file output of mHM
mhm.model.disable_output()
mhm.run.prepare()
# single domain run
mhm.run.prepare_domain()
while not mhm.run.finished():
    mhm.run.do_time_step()
    stream = mhm.get_variable("L11_QMOD")
    time = datetime(*mhm.run.current_time())
    if first_time:
        stream_acc = stream
        n = 1  # start averaging counter
        run_im = ax.imshow(stream, vmin=0, vmax=650)
        fig.colorbar(run_im, ax=ax, label="discharge in m³ / s")
        first_time = False
    else:
        # will reset for n=0
        stream_acc = (stream_acc * n + stream) / (n + 1)
        n += 1  # increase counter
    # plot each month
    if time.day == 1 and time.hour == 0:
        print("time:", time)
        n = 0  # reset counter
        run_im.set_data(stream_acc)
        fig.suptitle(f"{time}")
        fig.canvas.flush_events()
# finalize
mhm.run.finalize_domain()
mhm.run.finalize()
mhm.model.finalize()

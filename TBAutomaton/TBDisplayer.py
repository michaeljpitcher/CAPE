import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import gridspec


def create_movie(output_directory, grid_shape, movie_name, interval=400, legend=False, display=False):
    print "Collecting data..."
    contents_file = open(output_directory + "/contents.csv")
    spamreader = csv.reader(contents_file, delimiter=',')

    rows = []
    for row in spamreader:
        rows.append([float(r) for r in row])

    grids = []
    a = 0
    b = grid_shape[1]
    while b <= len(rows):
        grids.append(rows[a:b])
        a += grid_shape[1]
        b += grid_shape[1]

    time_steps = len(grids)

    fast_bacs = []
    fast_rest_bacs = []
    slow_bacs = []
    slow_rest_bacs = []
    rest_macs = []
    active_macs = []
    inf_macs = []
    chr_inf_macs = []
    t_cells = []
    caseum = []
    vessels = []


    for t in range(len(grids)):
        fast_bacs.append([])
        fast_rest_bacs.append([])
        slow_bacs.append([])
        slow_rest_bacs.append([])
        rest_macs.append([])
        active_macs.append([])
        inf_macs.append([])
        chr_inf_macs.append([])
        t_cells.append([])
        caseum.append([])
        grid = grids[t]
        for y in range(grid_shape[0]):
            for x in range(grid_shape[1]):

                if t == 0 and grid[x][y] == 1.5:
                    vessels.append((x,y))

                # FAST BAC
                if grid[x][y] == 1.0:
                    fast_bacs[t].append((x, y))
                # FAST BAC REST
                elif grid[x][y] == 1.25:
                    fast_rest_bacs[t].append((x, y))
                # SLOW BAC
                elif grid[x][y] == 2.0:
                    slow_bacs[t].append((x, y))
                # SLOW BAC REST
                elif grid[x][y] == 2.25:
                    slow_rest_bacs[t].append((x, y))
                # REST MAC
                elif grid[x][y] == 4.0:
                    rest_macs[t].append((x, y))
                # ACTIVE MAC
                elif grid[x][y] == 5.0:
                    active_macs[t].append((x, y))
                # INF MAC
                elif grid[x][y] == 6.0:
                    inf_macs[t].append((x, y))
                # CHR INF MAC
                elif grid[x][y] == 7.0:
                    chr_inf_macs[t].append((x, y))
                # T CELL
                elif grid[x][y] == 3.0:
                    t_cells[t].append((x, y))
                # CASEUM
                elif grid[x][y] == 100.0:
                    caseum[t].append((x, y))

    def update_plot(time_step):
        plt.clf()
        if legend:
            gs = gridspec.GridSpec(2, 1, height_ratios=[6, 1])
            plt.subplot(gs[0])

        plt.axis([0, grid_shape[0], grid_shape[1], 0])
        plt.xticks([])
        plt.yticks([])
        plt.suptitle("TB Automaton", fontsize=14, fontweight='bold')
        plt.title('Time = ' + str(time_step) + " hours", fontsize=10)

        bv = plt.scatter([v[1] for v in vessels], [v[0] for v in vessels],
                         s=20, color='red', marker="D")  # RED Blood Vessels
        fb = plt.scatter([fb[1] for fb in fast_bacs[time_step]], [fb[0] for fb in fast_bacs[time_step]],
                         s=1, color='#0F63AE')  # BLUE Fast Bacteria
        frb = plt.scatter([fbr[1] for fbr in fast_rest_bacs[time_step]], [fbr[0] for fbr in fast_rest_bacs[time_step]],
                          s=4, color='#20437c', marker="D")  # DEEP BLUE Fast resting Bacteria
        sb = plt.scatter([sb[1] for sb in slow_bacs[time_step]], [sb[0] for sb in slow_bacs[time_step]],
                         s=1, color='#851f98')  # PURPLE Slow Bacteria
        srb = plt.scatter([sbr[1] for sbr in slow_rest_bacs[time_step]], [sbr[0] for sbr in slow_rest_bacs[time_step]],
                          s=4, color='#490746', marker="D")  # DEEP PURPLE Slow Resting Bacteria
        rm = plt.scatter([rm[1] for rm in rest_macs[time_step]], [rm[0] for rm in rest_macs[time_step]],
                         color='#168964', marker=(5, 1))  # GREEN Resting macrophages
        am = plt.scatter([am[1] for am in active_macs[time_step]], [am[0] for am in active_macs[time_step]],
                         color='#00ff45', marker=(5, 1))  # BRIGHT GREEN Active Macrophages
        im = plt.scatter([im[1] for im in inf_macs[time_step]], [im[0] for im in inf_macs[time_step]],
                         color='#F1BC41', marker=(5, 1))  # GOLD Infected Macrophages
        cim = plt.scatter([cim[1] for cim in chr_inf_macs[time_step]], [cim[0] for cim in chr_inf_macs[time_step]],
                          color='#77643a', marker=(5, 1))  # BROWN Chronically Infected Macrophages
        tc = plt.scatter([tc[1] for tc in t_cells[time_step]], [tc[0] for tc in t_cells[time_step]],
                         color='#f9c7ed')  # PINK T-cells
        ca = plt.scatter([c[1] for c in caseum[time_step]], [c[0] for c in caseum[time_step]],
                         color='#000000')  # BLACK Caseum

        if legend:
            plt.subplot(gs[1])
            plt.xticks([])
            plt.yticks([])
            plt.legend((bv, fb, frb, sb, srb, rm, am, im, cim, tc, ca),
                       ("Blood vessel", "Fast bacterium", "Fast resting bacterium",
                        "Slow bacterium", "Slow resting bacterium", "Resting macrophage", "Active macrophage",
                        "Infected macrophage", "Chr. Infected macrophage", "T-cell", "Caseum"), scatterpoints=1,
                       loc='center', ncol=3, fontsize=8)

    # DISPLAY
    print "Creating animation..."
    fig = plt.figure()
    ani = animation.FuncAnimation(fig, update_plot, frames=xrange(time_steps), interval=interval, blit=False)
    ani.save(output_directory + "/" + movie_name + ".mp4", writer='ffmpeg_file')
    if display:
        plt.show()

if __name__ == '__main__':
    create_movie("../../../../Comparison/DATA/IV1/MINE2/3", (101,101), "TB")
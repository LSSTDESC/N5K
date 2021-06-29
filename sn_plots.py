import matplotlib.pyplot as plt

import numpy as np

if __name__ == "__main__":
    names = ["CCL", "Levin"]

    for name in names:
        comp_vs_bm = np.load(f"tests/tester_comp_{name}.npz")

        ell = comp_vs_bm["ls"]
        plt.semilogx(ell, comp_vs_bm["sn_per_l"], label=f"SN {name}")
        plt.semilogx(ell, np.sqrt(np.cumsum(comp_vs_bm["sn_per_l"]**2)-comp_vs_bm["sn_per_l"][0]**2), label=f"Cumulative SN {name}")

    plt.axhline(0, c="k", lw=1)
    
    plt.legend(frameon=False)
    plt.xlabel("$\\ell$")
    plt.ylabel("SN")
    plt.savefig("plots/sn_vs_ell.png")

    plt.figure()

    comp_vs_bm = np.load(f"tests/tester_comp_CCL.npz")
    ell = comp_vs_bm["ls"]
    err = np.concatenate([comp_vs_bm[f"cl_{t}_err"] for t in ["gg", "gs", "ss"]], axis=0)
    cl = np.concatenate([comp_vs_bm[f"cl_{t}"] for t in ["gg", "gs", "ss"]], axis=0)
    cl_bm = np.concatenate([comp_vs_bm[f"cl_{t}_bm"] for t in ["gg", "gs", "ss"]], axis=0)

    plt.imshow((cl-cl_bm)/err,
               extent=(ell[0], ell[-1], 0, cl.shape[0]),
               origin="lower")
    plt.colorbar(label="$(C_\\ell - C^{\\rm BM}_\\ell)/\\sigma$")

    plt.axhline(comp_vs_bm[f"cl_gg"].shape[0], c="k", lw=1)
    plt.axhline(comp_vs_bm[f"cl_gg"].shape[0]+comp_vs_bm[f"cl_gs"].shape[0], c="k", lw=1)

    plt.xscale("log")
    plt.xlabel("$\\ell$")

    plt.savefig("plots/err_vs_ell_vs_bin.png")


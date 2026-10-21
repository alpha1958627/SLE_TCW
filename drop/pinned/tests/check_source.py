"""Guard explicit y stencils; do not rewrite rotated foreach_face stencils."""
from pathlib import Path


def check(source):
    text = Path(source).read_text(encoding="utf-8")
    block = text.split("static void relax_hp001_coupled", 1)[1].split("static void hp001_fail", 1)[0]
    assert "M.y[1]" not in block
    assert block.count("M.y[0,1]") == 2
    noise = text.split("void apply_split_stochastic_flux", 1)[1]
    yface = noise.split("foreach_face(y)", 1)[1].split("boundary", 1)[0]
    assert "f_old[0,-1]" in yface and "f_old[-1]" not in yface
    volume = text.split("void update_volume_stats", 1)[1].split("int hp001_height_invalid", 1)[0]
    assert "MPI_Allreduce" not in volume
    assert "coeff_varphi != 0." in text
    return True


if __name__ == "__main__":
    import sys
    check(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parents[1]/"contact_angle.c")
    print("source stencil guards PASS")

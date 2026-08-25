"""Tests of the hadronic-final-state layer (polligen.hfs)."""

import numpy as np
import pytest

from polligen import hfs, recopseudo as rp
from polligen.sample import InclusiveSampler
from polligen.xsec import InclusiveKernel
from polli_fastsim import beams, fom
from polli_fastsim.polarized import toy_b1


@pytest.fixture(scope="module")
def toy_sample():
    rng = np.random.default_rng(7)
    cfg = beams.default_configs("6Li")[1]
    toy = hfs.ToyHFS(cfg.electron_energy, cfg.ion_momentum_per_nucleon)
    n = 4000
    x = 10 ** rng.uniform(-2.5, -0.7, n)
    q2 = 10 ** rng.uniform(0.0, 1.5, n)
    y = q2 / (toy.s * x)
    ok = (y > 0.005) & (y < 0.95)
    x, q2 = x[ok], q2[ok]
    smp = toy.generate(x, q2, rng.uniform(0, 2 * np.pi, x.size), rng)
    return cfg, smp


def _theta_e(smp):
    return np.arccos(smp.kp[:, 3] / np.sqrt((smp.kp[:, 1:] ** 2).sum(axis=1)))


def test_toy_conserves_four_momentum_and_truth_sums(toy_sample):
    cfg, smp = toy_sample
    ev = smp.event_index()
    tot = np.stack([np.bincount(ev, smp.p4[:, c], smp.n_events) for c in range(4)], axis=1)
    k = np.array([cfg.electron_energy, 0.0, 0.0, -cfg.electron_energy])
    P = np.array([cfg.ion_momentum_per_nucleon, 0.0, 0.0,
                  cfg.ion_momentum_per_nucleon])   # massless target (ToyHFS default)
    resid = tot + smp.kp - (k + P)[None, :]
    assert np.abs(resid).max() < 1e-6
    d_sig, d_pt = hfs.truth_kinematics_check(smp)
    # Sigma_true = 2 E_e y holds up to the target-mass term (E_N - p_N)
    assert d_sig < 2e-3 and d_pt < 1e-6
    assert np.all(np.diff(smp.offsets) >= 2)
    assert (smp.pid == 22).mean() > 0.1          # pi0 -> gamma gamma present


def test_kinematic_methods_are_exact_with_a_perfect_response(toy_sample):
    cfg, smp = toy_sample
    r = hfs.HadronResponse(perfect=True).reconstruct_sums(smp)
    kin = hfs.hadronic_kinematics(r["sigma"], r["ptx"], r["pty"], smp.kp[:, 0],
                                  _theta_e(smp), cfg.electron_energy, smp.s)
    for key in ("y_jb", "y_sigma", "y_da"):
        assert np.nanmax(np.abs(kin[key] / smp.y - 1.0)) < 5e-3, key
    assert np.nanmax(np.abs(kin["x_mixed"] / smp.x - 1.0)) < 5e-3
    assert np.nanmax(np.abs(kin["q2_e"] / smp.q2 - 1.0)) < 1e-9


def test_response_loses_energy_outside_acceptance_only(toy_sample):
    cfg, smp = toy_sample
    rng = np.random.default_rng(3)
    loose = hfs.HadronResponse(eta_track=50.0, eta_cal=50.0, pt_min_track=0.0,
                               pt_turn=1e-6, eff_track=1.0, e_min_photon=0.0,
                               e_min_nhad=0.0, eff_nhad=1.0, noise_sigma=0.0)
    r = loose.reconstruct_sums(smp, rng=rng)
    f = r["sigma"] / r["sigma_true"]
    # everything detected: only Gaussian smearing remains, centred on 1
    assert abs(np.median(f) - 1.0) < 0.02
    tight = hfs.HadronResponse(noise_sigma=0.0)
    r2 = tight.reconstruct_sums(smp, rng=rng)
    f2 = r2["sigma"] / r2["sigma_true"]
    assert np.median(f2) < np.median(f)          # thresholds and coverage lose Sigma
    assert np.percentile(f2, 99) < 1.5          # smearing tails only


def test_noise_is_additive_and_zero_mean(toy_sample):
    cfg, smp = toy_sample
    rng = np.random.default_rng(5)
    a = hfs.HadronResponse(noise_sigma=0.0).reconstruct_sums(smp, rng=np.random.default_rng(5))
    b = hfs.HadronResponse(noise_sigma=0.05).reconstruct_sums(smp, rng=np.random.default_rng(5))
    d = b["sigma"] - a["sigma"]
    assert abs(d.mean()) < 0.005 and abs(d.std() - 0.05) < 0.005


def test_library_transfer_reproduces_the_response_statistics(toy_sample):
    cfg, smp = toy_sample
    rng = np.random.default_rng(11)
    resp = hfs.HadronResponse(noise_sigma=0.0)
    lib = hfs.HFSLibrary(smp, resp, nx=12, nq2=8, rng=rng)
    frac, med = lib.coverage()
    assert frac > 0.3 and med >= 1
    f, r, dphi = lib.transfer(smp.x, smp.q2, rng)
    assert f.shape == smp.x.shape
    assert abs(np.median(f) - np.median(lib.f_sigma)) < 0.05
    # empty cells are redirected to populated ones
    f2, _, _ = lib.transfer(np.array([1e-4, 0.9]), np.array([1.0, 30.0]), rng)
    assert np.all(np.isfinite(f2))


def test_hfs_response_hooks_into_recopseudo():
    rng = np.random.default_rng(2)
    cfg = beams.default_configs("6Li")[1]
    analysis = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    sampler = InclusiveSampler(kern, cfg, rp.generator_scenario(analysis),
                               nx=20, nq2=15, q2_range=(0.7, 2e3))
    smp = hfs.toy_library_sample(sampler, 6000, rng)
    lib = hfs.HFSLibrary(smp, hfs.HadronResponse(), nx=20, nq2=15, rng=rng)
    model = rp.RecoModel(y_source="hfs", hadronic_method="sigma")
    resp = rp.RecoResponse(sampler, model, n_mc_per_cell=20, rng=rng,
                           hfs=hfs.HFSResponse(lib, method="sigma"))
    assert resp.y_reco.shape == resp.y.shape
    ok = resp.eff > 0
    ratio = resp.y_reco[ok] / resp.y[ok]
    # the hadronic y is centred near the truth and has a finite spread
    assert 0.5 < np.median(ratio) < 1.2
    assert 0.02 < np.std(np.log(ratio)) < 1.5
    assert resp.sigma_capture is not None and resp.sigma_capture.shape == resp.y.shape
    with pytest.raises(ValueError):
        rp.RecoResponse(sampler, rp.RecoModel(y_source="hfs"), n_mc_per_cell=5, rng=rng)


def test_sample_round_trip(tmp_path, toy_sample):
    cfg, smp = toy_sample
    p = tmp_path / "s.npz"
    smp.save(p)
    back = hfs.HFSSample.load(p)
    assert back.n_events == smp.n_events and back.pid.size == smp.pid.size
    assert np.allclose(back.p4, smp.p4) and back.meta.get("generator") == "ToyHFS"
    cat = hfs.HFSSample.concatenate([smp, back])
    assert cat.n_events == 2 * smp.n_events and cat.offsets[-1] == 2 * smp.offsets[-1]

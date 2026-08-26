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


# --- energy-scale nuisances and the beam guard (plans/08 A5) ----------------

def _library(cfg, smp, noise=0.05):
    hr = hfs.HadronResponse(noise_sigma=noise)
    return hfs.HFSLibrary(smp, hr, nx=24, nq2=18,
                          rng=np.random.default_rng(5))


def test_hadronic_scale_is_applied_after_the_noise(toy_sample):
    """The calibration multiplier belongs to HFSResponse, not to
    HadronResponse: the library's captured fraction f_sigma is built with
    the response's own parameters (noise off), so a scale there would be
    baked into f and applied twice."""
    cfg, smp = toy_sample
    lib = _library(cfg, smp)
    x, q2 = smp.x[:200], smp.q2[:200]
    s_nn = 4.0 * cfg.electron_energy * cfg.ion_momentum_per_nucleon
    y = np.clip(q2 / (s_nn * x), 1e-3, 0.9)
    e_p = np.full(x.size, 9.0)
    th = np.full(x.size, 3.0)
    args = (x, q2, y, e_p, th, cfg.electron_energy, s_nn)
    a = hfs.HFSResponse(lib).hadronic(*args, np.random.default_rng(1))
    b = hfs.HFSResponse(lib, scale=1.05).hadronic(*args,
                                                  np.random.default_rng(1))
    np.testing.assert_allclose(b["sigma"], 1.05 * a["sigma"], rtol=1e-12)
    np.testing.assert_allclose(a["f_sigma"], b["f_sigma"], rtol=1e-12)
    assert hfs.HFSResponse(lib, scale=1.0).scale == 1.0


def test_concatenate_refuses_to_merge_across_beam_energies(toy_sample):
    cfg, smp = toy_sample
    other = hfs.HFSSample(smp.offsets, smp.pid, smp.charge, smp.p4, smp.x,
                          smp.q2, smp.y, smp.kp, smp.weight,
                          smp.e_energy * 1.8, smp.p_per_nucleon,
                          {"target": "n"})
    with pytest.raises(ValueError, match="beam energies"):
        hfs.HFSSample.concatenate([smp, other])
    same = hfs.HFSSample(smp.offsets, smp.pid, smp.charge, smp.p4, smp.x,
                         smp.q2, smp.y, smp.kp, smp.weight, smp.e_energy,
                         smp.p_per_nucleon, {"target": "n"})
    merged = hfs.HFSSample.concatenate([smp, same])
    assert merged.n_events == 2 * smp.n_events
    assert len(merged.meta["merged"]) == 2          # p+n merge is allowed


def test_reco_response_refuses_a_library_from_other_beams(toy_sample):
    """The coupling point, so the mistake is caught on any call path."""
    cfg, smp = toy_sample
    lib = _library(cfg, smp)
    resp = hfs.HFSResponse(lib)
    resp.check_beams(cfg.electron_energy, cfg.ion_momentum_per_nucleon)
    with pytest.raises(ValueError, match="not transferable"):
        resp.check_beams(cfg.electron_energy, 137.5)
    other_cfg = beams.default_configs("6Li")[2]
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    scen = rp.generator_scenario(fom.Scenario())
    sampler = InclusiveSampler(kern, other_cfg, scen, nx=8, nq2=6,
                               q2_range=(0.7, 2e3))
    with pytest.raises(ValueError, match="not transferable"):
        rp.RecoResponse(sampler, rp.RecoModel(y_source="hfs"),
                        n_mc_per_cell=5, rng=np.random.default_rng(1),
                        hfs=resp)


def test_energy_scale_levers_of_the_sigma_method(toy_sample):
    """d ln x / d ln E' = 2 - y and d ln x / d ln (hadronic scale) =
    -(1 - y): the ELECTRON scale is the bigger lever, about twice the
    hadronic one."""
    cfg, smp = toy_sample
    lib = _library(cfg, smp)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    scen = fom.Scenario()
    gen = rp.generator_scenario(scen)

    def build(e_scale=1.0, h_scale=1.0):
        sampler = InclusiveSampler(kern, cfg, gen, nx=20, nq2=15,
                                   q2_range=(0.7, 2e3))
        model = rp.RecoModel(q2_min=scen.q2_min, y_min=scen.y_min,
                             y_max=scen.y_max, w2_min=scen.w2_min,
                             eta_min=scen.eta_min, eta_max=scen.eta_max,
                             e_prime_min=scen.e_prime_min, y_source="hfs",
                             e_scale=e_scale)
        return rp.RecoResponse(sampler, model, n_mc_per_cell=60,
                               rng=np.random.default_rng(7),
                               hfs=hfs.HFSResponse(lib, scale=h_scale))

    base = build()
    for kw, expected in ((dict(e_scale=1.01), lambda y: 2.0 - y),
                         (dict(h_scale=1.01), lambda y: -(1.0 - y))):
        alt = build(**kw)
        ok = (np.isfinite(base.x_reco) & np.isfinite(alt.x_reco)
              & (base.x_reco > 0) & (alt.x_reco > 0))
        lever = np.median(np.log(alt.x_reco[ok] / base.x_reco[ok])
                          / np.log(1.01))
        assert lever == pytest.approx(expected(np.median(base.y[ok])),
                                      abs=0.05)


def test_parametrized_hadronic_y_understates_the_electron_lever(toy_sample):
    """A limitation of the Gaussian stand-in worth stating: it smears the
    TRUE y and never sees E', so the electron energy scale moves x only
    through Q2 (lever 1) instead of 2 - y."""
    cfg, _ = toy_sample
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    scen = fom.Scenario()
    gen = rp.generator_scenario(scen)

    def build(e_scale=1.0):
        sampler = InclusiveSampler(kern, cfg, gen, nx=20, nq2=15,
                                   q2_range=(0.7, 2e3))
        model = rp.RecoModel(q2_min=scen.q2_min, y_min=scen.y_min,
                             y_max=scen.y_max, w2_min=scen.w2_min,
                             eta_min=scen.eta_min, eta_max=scen.eta_max,
                             e_prime_min=scen.e_prime_min, e_scale=e_scale)
        return rp.RecoResponse(sampler, model, n_mc_per_cell=60,
                               rng=np.random.default_rng(7))

    base, alt = build(), build(1.01)
    ok = (np.isfinite(base.x_reco) & np.isfinite(alt.x_reco)
          & (base.x_reco > 0) & (alt.x_reco > 0))
    lever = np.median(np.log(alt.x_reco[ok] / base.x_reco[ok]) / np.log(1.01))
    assert lever == pytest.approx(1.0, abs=1e-6)

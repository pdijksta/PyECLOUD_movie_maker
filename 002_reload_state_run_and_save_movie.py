from PyECLOUD.buildup_simulation import BuildupSimulation


sim = BuildupSimulation()
sim.load_state('simulation_state_0.pkl')
sim.pyeclsaver.flag_video = True
sim.pyeclsaver.flag_sc_video = True
sim.run(t_end_sim=sim.beamtim.tt_curr+50e-9)

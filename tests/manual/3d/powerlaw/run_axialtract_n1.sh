#!/bin/bash

pylith axialtraction_powerlaw_n1_norefstate.cfg axialtraction_powerlaw_n1_norefstate_dt0.01_onecell_hex.cfg > axialtraction_powerlaw_n1_norefstate_dt0.01_onecell_hex.log 2>&1
pylith axialtraction_powerlaw_n1_norefstate.cfg axialtraction_powerlaw_n1_norefstate_dt0.02_onecell_hex.cfg > axialtraction_powerlaw_n1_norefstate_dt0.02_onecell_hex.log 2>&1
pylith axialtraction_powerlaw_n1_norefstate.cfg axialtraction_powerlaw_n1_norefstate_dt0.05_onecell_hex.cfg > axialtraction_powerlaw_n1_norefstate_dt0.05_onecell_hex.log 2>&1
pylith axialtraction_powerlaw_n1_norefstate.cfg axialtraction_powerlaw_n1_norefstate_dt0.1_onecell_hex.cfg > axialtraction_powerlaw_n1_norefstate_dt0.1_onecell_hex.log 2>&1
pylith axialtraction_powerlaw_n1_norefstate.cfg axialtraction_powerlaw_n1_norefstate_dt0.2_onecell_hex.cfg > axialtraction_powerlaw_n1_norefstate_dt0.2_onecell_hex.log 2>&1

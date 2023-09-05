#!/usr/bin/python

import scrublet as scr

def score_doublets(mtx, doublet_rate=0.1):
    scrub = scr.Scrublet(mtx.T, expected_doublet_rate=doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    return doublet_scores

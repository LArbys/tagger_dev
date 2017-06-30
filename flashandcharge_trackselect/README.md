# Improving the Tagger by using Flash and Charge selection metrics

One major challenge of the tagger is determining when the thrumu or stopmu tracks are "bad".

Bad tracks include

* ThruMu tracks that only cover a portion of the true track, due to false
  end points located in the middle (and not the ends) of tracks
* Spurious ThruMu tracks from the A* algorithm. Usually caused by passing
  through areas with a lot of dead regions.  Causes the 3D path finder to project onto regions
  on each of the planes over different tracks
* StopMu tracks that start in the middle of tracks, again, due to the spurious end point finding


This development repository aims to include flash matching and calorimetry much earlier to
mitigate these issues.

* Use flash-matching on through-going muons early, in order to reject, short, early-stopping
  track candidates
* Match thrumu and stop mu tracks with anode/cathode crossing flash ends. Finding the proper
  match early allows us to reject spurious anode/cathode candidate end points early
* Compare the charge profile and totals between the 3 planes to reject spurious tracks.
  This depends on the quality of the calorimetry modeling between planes.  This might not
  work as the current deconvolution might not be stable/calibrated well enough between planes
  and between data-to-MC.

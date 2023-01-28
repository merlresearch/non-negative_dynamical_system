<!--
Copyright (C) 2012,2013,2023 Mitsubishi Electric Research Laboratories (MERL)

SPDX-License-Identifier: AGPL-3.0-or-later
-->
# Non-negative Dynamical System model (NDS)

## Summary

MATLAB source code for the ICASSP 2013 paper **Non-negative Dynamical System with Application to Speech and Audio** by Cédric Févotte, Jonathan Le Roux, and John R. Hershey.

[Please click here to read this paper.](https://www.merl.com/publications/TR2013-021)

## Usage

The code should run without requiring installation.

We provide a [demo script](nds_demo.m) and associated audio samples to let users test the software. The NDS model is trained on around 13 s of speech by a female speaker of the Librispeech train-clean-360 dataset (`training_speech.wav`), and tested on the mixture (`noisy.wav`) of another female speaker of Librispeech train-clean-360 with a helicopter sound.

## Citation

If you use any part of this code for your work, we ask that you include the following citation:
```
@inproceedings{Fevotte2013ICASSP05,
    author = {F\'{e}votte, C\'{e}dric and {Le Roux}, Jonathan and Hershey, John R.},
    title = {Non-negative Dynamical System with Application to Speech and Audio},
    booktitle = {Proc. IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP)},
    year = 2013,
    month = may,
    url = {https://www.merl.com/publications/TR2013-021}
}
```

## Related Links

Please refer to the [Non-Negative Dynamical System With Application to Speech And Audio demo page](https://www.merl.com/demos/speech-enhancement) for demo videos and audio samples.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for our policy on contributions.

## License

Released under `AGPL-3.0-or-later` license, as found in the [LICENSE.md](LICENSE.md) file.

All files, except as listed below:

```
Copyright (c) 2012,2013,2023 Mitsubishi Electric Research Laboratories (MERL).

SPDX-License-Identifier: AGPL-3.0-or-later
```

Files in `auxfun/` are licensed under `GPL-3.0-or-later` (see [LICENSES/GPL-3.0-or-later.md](LICENSES/GPL-3.0-or-later.md)).
```
Copyright Cedric Fevotte (CNRS), 2012.

SPDX-License-Identifier: GPL-3.0-or-later
```

Sample file `training_speech.wav` is licensed under Creative Commons Attribution 4.0 International License `CC-BY-4.0` (see [LICENSES/CC-BY-4.0.txt](LICENSES/CC-BY-4.0.txt)), as it is derived from the Librispeech corpus:

```
Copyright (c) 2014 by Vassil Panayotov

SPDX-License-Identifier: CC-BY-4.0
```

Sample file `noisy.wav` is derived from [Sound Textures](https://mcdermottlab.mit.edu/McDermott_Simoncelli_2011_168_Sound_Textures.zip), introduced in "Sound Texture Perception via Statistics
of the Auditory Periphery: Evidence from Sound Synthesis", Josh H. McDermott and Eero P. Simoncelli, Neuron, 71(5), pp.926-940 and also from the Librispeech corpus:

```
Copyright (c) 2023 Mitsubishi Electric Research Laboratories (MERL).
Copyright (c) 2014 by Vassil Panayotov
Copyright (c) 2011 Josh H. McDermott and Eero P. Simoncelli

SPDX-License-Identifier: AGPL-3.0-or-later
```

# MacRes - short for Mac Resampler

Fork of https://github.com/ohac/tn_fnds.

This is an executable designed to be called by https://github.com/titinko/utsu (UTSU).
Despite the name, it can be compiled for Windows, Mac, or Linux using the Makefile.

The Windows/Linux versions of this resampler can also be used with the UTAU program:
https://en.wikipedia.org/wiki/Utau.


# Compilation

To build macres just do
```
cd src
make
```


# Use

```
macres inputfile outputfile notenum [velocity flags offset_ms notelength fixedlength end intensity modulation tempo pitchbends]
```

* Inputfile: Wav input file
* Outputfile: Wav output file
* notenum: Pitch to resample. Could use C0-B7 . Use C# or Db for flat and bemol.
* velocity: Consonant velocity, in percent. 100 uses the original length, 50 half time, 200 double velocity
* Flags: Use ? for no flags
  * B0 - B100 : Breath . Default=50
  * b0 - b100 : Consonant strength. Default=0
  * t<XXX> : transpose XXX/120 semitones. For example, t120 traspose the note +1oct. 
  * g-100 - g100 : Gender flag
  * W-1 unvoiced W0 no effect W50-5000 frequency enforcement
  * A0 - A100: Correct volume of combined
  * O0 - O100: Voice strength. Change the voice "brightness". Specifying + values suppresses lower frequencies and amplifies higher one
  * e : Change the vowel stretching method to UTSU
* Offset: Offset to start the note in ms 
* Note length: Total note duration in ms
* Fixed Length: Consonant section in ms. Not used on loop when we extend the note 
* End: Offset from the end of the file to remove.
* Intensity: Volume of the note. 100 is same as original
* Modulation: Variations around the f0. 0=no modulation, perfect tone (autotune-like). 100= original modulations around the average f0. 200= double of modulation around the average.
* Tempo: Tempo of the song, on BPM. By default, 120. Should be prepended by 'T'
* Pitchbend: String that represents changes in pitch. See *Pitchbend specification* for more information

# Pitchbend specification

You can change the pitch of the note using a pitchbend string. On every beat you have 96 pitch step, so for *s* seconds on $bpm$ tempo, you have $\frac{96 bmp*s}{60}$ steps available.
Every pitch change is declared using a base-64 notation, based on couples of letters, first letter more significant. The availables values are *A-Za-z0-9+/* .Minimum is 'A'=0, 'Z'=25 an so on, until '/'=63. 

| s | v    |
|---|------|
| A |    0 |
| B |    1 |
| C |    2 |
| D |    3 |
|   |  ... |
| Z |   25 |
| a |   26 |
| b |   27 |
| c |   28 |
| d |   29 |
| e |   30 |
| f |   31 |
| g |   32 |
| h |   33 |
| i |   34 |
| j |   35 |
| k |   36 |
| . |   .. |
| z |   .. |
| 4 |   56 |
| 5 |   57 |
| 6 |   58 |
| 7 |   59 |
| 8 |   60 |
| 9 |   61 |
| + |   62 |
| / |   63 |
------------

A change of 1 semitone is represented by a difference of 100, so is possible to pitch about 20 semitones (octave and a fifth). Values from 0  (*AA*) upto 2047 (*e/*) are considered positive. From 2048(*fA*) upto 4196 (*/*) are considered negatives, using $x-4096$. 

So:


* No change: AA (0)
* 1 semitone up: Bk (100)
* 2 semitones up: DI (200)
* 1 semitone down: +c (3996-4096= -100)
* 2 semitones down: 84 (3896-4096 = -200)

To repeat a value, you use the format #TIMES# after the note. So, to maintain the base note for  100 steps, you could use AA#100#


# Examples

Use the file 'mam.wav' in *assets* dir as example. Is mainly a C3 (130.81 Hz), with some problems of tuning at the first and last 'm'. The useful part of the note is between 750 ms and 2500ms.
Note the differences on amplitude of consonants and vowels, created to help to understand the changes of the resampler.

![Mam spectrograph](assets/mam.png)

A fairly vanilla setting, with a 5000ms note, 100% velocity of consonant, 100% volume, 100% of original modulation and taking off the noise will be:
```
macres assets/mam.wav test.wav C3 100 ? 750 2000 200 1250 100 100
```

![Eg.1](assets/eg_1.png)

Note the 200ms for the consonant. If we remove the consonant, the consonant will be played again when original sample length is shorter than required note:

```
macres assets/mam.wav test.wav C3 100 ? 750 5000 000 1250 100 100
```

![Eg.2](assets/eg_2.png)

To generate a perfect pitch, we could set modulation to 0:
```
macres assets/mam.wav test.wav C3 100 ? 750 5000 200 1250 100 0
```
![Eg.3](assets/eg_3.png)


Using a tempo of 60, we have 96 pitch per second. On this 2 second note, we have 192 steps and we try a little up and down glissando
```
macres assets/mam.wav test.wav C3 100 ? 750 2000 200 1250 100 0 T60 AA#26#Bk#20#DI#20#Bk#20#AA#20#+c#20#84#20#+c#20#AA#26#
```
 
![Eg.4](assets/eg_4.png)

# Implementation

This resampler uses the  [World](https://github.com/mmorise/World) library to implements the synthesis.
Specifically, uses some (relatively old) specs:

* DIO: estimate F0 contour
* Platinum: Method to extract excitation signals for voice synthesis system based on PLATINUM


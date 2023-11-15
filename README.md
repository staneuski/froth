# froth
`froth` is OpenFOAM addons (including libraries) and `$WM_PROJECT_SITE` case templates.

## addons
Add custom #includeEtc dictionaries by setting `WM_PROJECT_SITE`
```sh
export WM_PROJECT_SITE=$HOME/Documents/froth
```

## ParaView
Animate ParaView state using `pvbatch`:
```sh
field=U
pvbatch $WM_PROJECT_SITE/etc/visualise.py state.pvsm animate /tmp/frames/$field.png --dict $WM_PROJECT_SITE/etc/caseDicts/postProcessing/visualisation/animation.png.json &&
ffmpeg -y -framerate 10 -pattern_type glob -i "$c$s$t/postProcessing/frames/$field.*.png" -qscale:v 8 -codec:a libvorbis postProcessing/U.ogv &&
rm -rf postProcessing/frames/$field.*.png
```

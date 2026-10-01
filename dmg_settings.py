"""dmgbuild settings for PyChem. Invoked with `-D app=dist/PyChem.app`."""

from pathlib import Path

try:
    APP = Path(app)  # noqa: F821 — injected by dmgbuild via -D app=...
except NameError:
    APP = Path("dist/PyChem.app")

format = "UDZO"
size = None
files = [str(APP)]
symlinks = {"Applications": "/Applications"}
icon = "assets/icon.icns"
badge_icon = icon
icon_locations = {
    APP.name: (140, 170),
    "Applications": (420, 170),
}
background = None
window_rect = ((200, 200), (580, 360))
default_view = "icon-view"
show_icon_preview = False
include_icon_view_settings = True
icon_view_settings = {
    "text_size": 14,
    "icon_size": 96,
    "arrange_by": None,
}

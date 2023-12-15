from bs4 import BeautifulSoup
import re
import csv

keypoints_vals = {
    "xmax" : 520.0,
    "ymax" : 0.124,
    "xmin" : 0.0,
    "ymin" : 0.0
}


with open('./literature/gas_sorption_data.svg', 'r') as f:
    data = f.read()


def getStartPoint(path: str) -> (float, float):
    m = re.search(r"(-*\d*.\d*),(-*\d*.\d*)", path)
    xpos, ypos = float(m.group(1)), float(m.group(2))
    return xpos, ypos
 
def getBoundryOfRectPath(path: str) -> (float, float, float, float):
    xpos, ypos = getStartPoint(path)
    dx, dy = 0, 0

    m = re.search(r"(H|h)\s*(-*\d*.\d*)", path)
    token = m.group(1)
    value = float(m.group(2))
    if token == "H":
        dx = value - xpos
    else:
        dx = value

    m = re.search(r"(V|v)\s*(-*\d*.\d*)", path)
    token = m.group(1)
    value = float(m.group(2))
    if token == "V":
        dy =value - ypos
    else:
        dy = value

    return xpos, ypos, dx, dy

svg_data = BeautifulSoup(data, "xml")
border_data = svg_data.find('path', {"inkscape:label": "border"})
print(border_data.get("d"))
print(getBoundryOfRectPath(border_data.get("d")))


keypoints_coord = {}

for key, value in keypoints_vals.items():
    path = svg_data.find('path', {"inkscape:label": key})
    if key.startswith("x"):
        keypoints_coord[key], _ = getStartPoint(path.get('d'))
    else:
        _, keypoints_coord[key] = getStartPoint(path.get('d'))

print(keypoints_coord)


def world2local(x, y) -> (float, float):
    x0, y0 = keypoints_coord["xmin"], keypoints_coord["ymin"]
    x1, y1 = keypoints_coord["xmax"], keypoints_coord["ymax"]

    dx = (keypoints_vals["xmax"] - keypoints_vals["xmin"]) / (x1 - x0)
    dy = (keypoints_vals["ymax"] - keypoints_vals["ymin"]) / (y1 - y0)
    
    result = (x - x0) * dx , (y - y0) * dy
    return result

points_data = svg_data.find('g', {"inkscape:label": "points_dV"})
points = points_data.findChildren("path", recursive=False)
datas = []
for point in points:
    path = point.get("d")
    boundry = getBoundryOfRectPath(path)
    x, y = world2local(
            boundry[0] + boundry[2] / 2.0, 
            boundry[1] + boundry[3] / 2.0
        )
    datas.append(
        [x, y]
    )

datas.sort(key=lambda e: e[0])
# print(datas)


with open('./data/pore_dv_data.csv', 'w') as csvfile:
    writer = csv.writer(csvfile,quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["pore width [A]", "V [cm^3/A/g]"])
    writer.writerows(datas)
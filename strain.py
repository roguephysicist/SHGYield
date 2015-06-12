import math

mode = 'normal'
DELTA = 10

if mode == 'big':
    Si1 = [0.00000000000000, 0.00000000000000, 15.99564782028059]
    Si2 = [4.18912260765629, 0.00000000000000, 14.51456931840277]
elif mode == 'normal':
    Si1 = [0.00000000000000, 0.00000000000000, 15.55132426971724]
    Si2 = [4.18912260765629, 0.00000000000000, 14.07024576783943]
elif mode == 'small':
    Si1 = [0.00000000000000, 0.00000000000000, 15.10700071915310]
    Si2 = [4.18912260765629, 0.00000000000000, 13.62592221727608]

H = 2.8013111110535096
L = 4.44323550563346

def distance(point1, point2):
    dist = math.sqrt(((point2[0] - point1[0]) ** 2) + (point2[2] - point1[1]) ** 2) + ((point2[2] - point1[2]) ** 2))
    return dist

percent = 1 + DELTA/100.0
hnew = distance(Si1, Si2) * percent
zdist = Si1[2] - Si2[2]
xdist = Si1[0] - Si2[0]
ang = math.atan(zdist/xdist) + math.pi
znew = hnew * math.sin(ang) + Si2[2]
xnew = hnew * math.cos(ang) + Si2[0]
print "{0:17.14f}   {1:17.14f}    {2:17.14f}".format(xnew, Si2[1], znew+H)
print "{0:17.14f}   {1:17.14f}    {2:17.14f}".format(xnew, Si2[1], znew)



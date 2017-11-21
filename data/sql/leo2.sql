Leo2 position in real units
RA:  168.3708
DEC: 22.1517

---------------------------
SELECT
 p.ra, p.dec, p.psfMag_g, psfMagErr_g, p.psfMag_r, p.psfMagErr_r, p.psfMag_i, p.psfMagErr_i
FROM PhotoObjAll AS p
WHERE
  p.ra BETWEEN 167.3708 AND 169.3708
  AND p.dec BETWEEN 21.1517 AND 23.1517
  AND dbo.fDistanceArcMinEq(168.2, 22.1517, p.ra, p.dec) < 30.
  AND p.clean = 1


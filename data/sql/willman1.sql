SELECT
 p.ra, p.dec, p.psfMag_g, psfMagErr_g, p.psfMag_r, p.psfMagErr_r, p.psfMag_i, p.psfMagErr_i
FROM PhotoObjAll AS p
WHERE
  p.ra BETWEEN 161.3375 AND 163.3375
  AND p.dec BETWEEN 50.05 AND 52.05
  AND dbo.fDistanceArcMinEq(162.3375, 51.05, p.ra, p.dec) < 30.
  AND p.clean=1

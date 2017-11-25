SELECT
 p.ra, p.dec, p.psfMag_g, psfMagErr_g, p.psfMag_r, p.psfMagErr_r, p.psfMag_i, p.psfMagErr_i
FROM PhotoObjAll AS p
WHERE
  p.ra BETWEEN 150.77 AND 152.77
  AND p.dec BETWEEN 15.08 AND 17.08
  AND dbo.fDistanceArcMinEq(151.77, 16.08, p.ra, p.dec) < 30.
  AND p.clean = 1


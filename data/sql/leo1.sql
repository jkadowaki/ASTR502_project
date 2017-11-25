SELECT
 p.ra, p.dec, p.psfMag_g, psfMagErr_g, p.psfMag_r, p.psfMagErr_r, p.psfMag_i, p.psfMagErr_i
FROM PhotoObjAll AS p
WHERE
  p.ra BETWEEN 151.12 AND 153.12
  AND p.dec BETWEEN 11.3 AND 13.3
  AND dbo.fDistanceArcMinEq(152.12, 12.3, p.ra, p.dec) < 30.
  AND p.clean = 1

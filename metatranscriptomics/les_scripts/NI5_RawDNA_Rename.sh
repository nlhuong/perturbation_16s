#!/bin/bash
# Les Dethlefsen 20 March 2018
# Add MIDs to start of No Intervention raw DNA read files
# Simplify full codes

# Note there are 4 files from DNA_Plate_10 that are not part of this work and should be dropped:
# rm AAC_CC_ST_Ellen_49_5d*
# rm AAC_CC_ST_Ellen_50_5d*
# rm AAC_CC_ST_Ellen_51_5d*
# rm AAC_CC_ST_Ellen_52_5d*

# Initial steps on 2 lines below necessary to reuse earlier renaming script
rename _NI_SW_ sw_ *_NI_SW_*
rename _NI_ST_ st_ *_NI_ST_*

# Following 6 lines needed due to inconsistent file names on DNA_Plate_12
rename CAC_ CACst_ CAC_??d*
rename CAD_ CADst_ CAD_??d*
rename CAD_ CADst_ CAD_?d*
rename CAM_ CAMst_ CAM_??d*
rename CAN_ CANst_ CAN_??d*
rename CAN_ CANst_ CAN_?d*

rename CAKsw_27d M1296_CAKsw_27d CAKsw_27d*
rename CAKsw_28d M1297_CAKsw_28d CAKsw_28d*
rename CAKsw_29d M1298_CAKsw_29d CAKsw_29d*
rename CAKsw_30d M1299_CAKsw_30d CAKsw_30d*
rename CAKsw_31d M1300_CAKsw_31d CAKsw_31d*
rename CAKsw_32d M1301_CAKsw_32d CAKsw_32d*
rename CAKsw_33d M1302_CAKsw_33d CAKsw_33d*
rename CAKsw_34d M1303_CAKsw_34d CAKsw_34d*
rename CAKsw_35d M1304_CAKsw_35d CAKsw_35d*
rename CAKsw_36d M1305_CAKsw_36d CAKsw_36d*
rename CAKsw_37d M1306_CAKsw_37d CAKsw_37d*
rename CAKsw_38d M1307_CAKsw_38d CAKsw_38d*
rename CAKsw_39d M1308_CAKsw_39d CAKsw_39d*
rename CAKsw_40d M1309_CAKsw_40d CAKsw_40d*
rename CAKsw_41d M1310_CAKsw_41d CAKsw_41d*
rename CAKsw_42d M1311_CAKsw_42d CAKsw_42d*
rename CAKsw_43d M1312_CAKsw_43d CAKsw_43d*
rename CAKsw_44d M1313_CAKsw_44d CAKsw_44d*
rename CAKsw_45d M1314_CAKsw_45d CAKsw_45d*
rename CAKsw_46d M1315_CAKsw_46d CAKsw_46d*
rename CAKsw_47d M1316_CAKsw_47d CAKsw_47d*
rename CAKsw_48d M1317_CAKsw_48d CAKsw_48d*
rename CAKsw_49d M1318_CAKsw_49d CAKsw_49d*
rename CAKsw_50d M1319_CAKsw_50d CAKsw_50d*
rename CAKsw_51d M1320_CAKsw_51d CAKsw_51d*
rename CAKsw_52d M1321_CAKsw_52d CAKsw_52d*
rename CAMsw_17d M1323_CAMsw_17d CAMsw_17d*
rename CAMsw_18d M1324_CAMsw_18d CAMsw_18d*
rename CAMsw_19d M1325_CAMsw_19d CAMsw_19d*
rename CAMsw_20d M1326_CAMsw_20d CAMsw_20d*
rename CAMsw_21d M1327_CAMsw_21d CAMsw_21d*
rename CAMsw_22d M1328_CAMsw_22d CAMsw_22d*
rename CAMsw_26d M1329_CAMsw_26d CAMsw_26d*
rename CAMsw_27d M1330_CAMsw_27d CAMsw_27d*
rename CAMsw_28d M1331_CAMsw_28d CAMsw_28d*
rename CAMsw_29d M1332_CAMsw_29d CAMsw_29d*
rename CAMsw_30d M1333_CAMsw_30d CAMsw_30d*
rename CAMsw_31d M1334_CAMsw_31d CAMsw_31d*
rename CAMsw_32d M1335_CAMsw_32d CAMsw_32d*
rename CAMsw_33d M1336_CAMsw_33d CAMsw_33d*
rename CAMsw_34d M1337_CAMsw_34d CAMsw_34d*
rename CAMsw_35d M1338_CAMsw_35d CAMsw_35d*
rename CAMsw_36d M1339_CAMsw_36d CAMsw_36d*
rename CAMsw_37d M1340_CAMsw_37d CAMsw_37d*
rename CAMsw_38d M1341_CAMsw_38d CAMsw_38d*
rename CAMsw_39d M1342_CAMsw_39d CAMsw_39d*
rename CAMsw_40d M1343_CAMsw_40d CAMsw_40d*
rename CAMsw_41d M1344_CAMsw_41d CAMsw_41d*
rename CAMsw_42d M1345_CAMsw_42d CAMsw_42d*
rename CAMsw_43d M1346_CAMsw_43d CAMsw_43d*
rename CAMsw_44d M1347_CAMsw_44d CAMsw_44d*
rename CAMsw_45d M1348_CAMsw_45d CAMsw_45d*
rename CAMsw_46d M1349_CAMsw_46d CAMsw_46d*
rename CAMsw_47d M1350_CAMsw_47d CAMsw_47d*
rename CAMsw_48d M1351_CAMsw_48d CAMsw_48d*
rename CAMsw_49d M1352_CAMsw_49d CAMsw_49d*
rename CAMsw_50d M1353_CAMsw_50d CAMsw_50d*
rename CAMsw_51d M1354_CAMsw_51d CAMsw_51d*
rename CANsw_16d M1356_CANsw_16d CANsw_16d*
rename CANsw_17d M1357_CANsw_17d CANsw_17d*
rename CANsw_18d M1358_CANsw_18d CANsw_18d*
rename CANsw_19d M1359_CANsw_19d CANsw_19d*
rename CANsw_20d M1360_CANsw_20d CANsw_20d*
rename CANsw_21d M1361_CANsw_21d CANsw_21d*
rename CANsw_22d M1362_CANsw_22d CANsw_22d*
rename CANsw_23d M1363_CANsw_23d CANsw_23d*
rename CANsw_24d M1364_CANsw_24d CANsw_24d*
rename CANsw_25d M1365_CANsw_25d CANsw_25d*
rename CANsw_26d M1366_CANsw_26d CANsw_26d*
rename CANsw_27d M1367_CANsw_27d CANsw_27d*
rename CANsw_28d M1368_CANsw_28d CANsw_28d*
rename CANsw_29d M1369_CANsw_29d CANsw_29d*
rename CANsw_30d M1370_CANsw_30d CANsw_30d*
rename CANsw_31d M1371_CANsw_31d CANsw_31d*
rename CANsw_32d M1372_CANsw_32d CANsw_32d*
rename CANsw_33d M1373_CANsw_33d CANsw_33d*
rename CANsw_34d M1374_CANsw_34d CANsw_34d*
rename CANsw_35d M1375_CANsw_35d CANsw_35d*
rename CANsw_36d M1376_CANsw_36d CANsw_36d*
rename CANsw_37d M1377_CANsw_37d CANsw_37d*
rename CANsw_38d M1378_CANsw_38d CANsw_38d*
rename CANsw_39d M1379_CANsw_39d CANsw_39d*
rename CANsw_40d M1380_CANsw_40d CANsw_40d*
rename CANsw_41d M1381_CANsw_41d CANsw_41d*
rename CANsw_42d M1382_CANsw_42d CANsw_42d*
rename CANsw_43d M1383_CANsw_43d CANsw_43d*
rename CANsw_44d M1384_CANsw_44d CANsw_44d*
rename CANsw_45d M1385_CANsw_45d CANsw_45d*
rename CANsw_46d M1386_CANsw_46d CANsw_46d*
rename CANsw_47d M1387_CANsw_47d CANsw_47d*
rename CANsw_48d M1388_CANsw_48d CANsw_48d*
rename CANsw_49d M1389_CANsw_49d CANsw_49d*
rename CANsw_50d M1390_CANsw_50d CANsw_50d*
rename PL9_B1d M1322_PL9_B1d PL9_B1d*
rename PL9_B2d M1355_PL9_B2d PL9_B2d*
rename PL9_B3d M1391_PL9_B3d PL9_B3d*
rename CAAst_10d M3417_CAAst_10d CAAst_10d*
rename CAAst_11d M3418_CAAst_11d CAAst_11d*
rename CAAst_12d M3419_CAAst_12d CAAst_12d*
rename CAAst_13d M3420_CAAst_13d CAAst_13d*
rename CAAst_14d M3421_CAAst_14d CAAst_14d*
rename CAAst_15d M3422_CAAst_15d CAAst_15d*
rename CAAst_16d M3423_CAAst_16d CAAst_16d*
rename CAAst_17d M3424_CAAst_17d CAAst_17d*
rename CAAst_18d M3425_CAAst_18d CAAst_18d*
rename CAAst_19d M3426_CAAst_19d CAAst_19d*
rename CAAst_1d M3408_CAAst_1d CAAst_1d*
rename CAAst_20d M3427_CAAst_20d CAAst_20d*
rename CAAst_21d M3428_CAAst_21d CAAst_21d*
rename CAAst_22d M3429_CAAst_22d CAAst_22d*
rename CAAst_23d M3430_CAAst_23d CAAst_23d*
rename CAAst_24d M3431_CAAst_24d CAAst_24d*
rename CAAst_25d M3432_CAAst_25d CAAst_25d*
rename CAAst_26_5d M3497_CAAst_26_5d CAAst_26_5d*
rename CAAst_26d M3433_CAAst_26d CAAst_26d*
rename CAAst_27_5d M3498_CAAst_27_5d CAAst_27_5d*
rename CAAst_27d M3434_CAAst_27d CAAst_27d*
rename CAAst_28_5d M3499_CAAst_28_5d CAAst_28_5d*
rename CAAst_28d M3435_CAAst_28d CAAst_28d*
rename CAAst_29d M3436_CAAst_29d CAAst_29d*
rename CAAst_2d M3409_CAAst_2d CAAst_2d*
rename CAAst_30d M3437_CAAst_30d CAAst_30d*
rename CAAst_31d M3438_CAAst_31d CAAst_31d*
rename CAAst_32d M3439_CAAst_32d CAAst_32d*
rename CAAst_33d M3440_CAAst_33d CAAst_33d*
rename CAAst_34d M3441_CAAst_34d CAAst_34d*
rename CAAst_35d M3442_CAAst_35d CAAst_35d*
rename CAAst_36d M3443_CAAst_36d CAAst_36d*
rename CAAst_37d M3444_CAAst_37d CAAst_37d*
rename CAAst_38d M3445_CAAst_38d CAAst_38d*
rename CAAst_39d M3446_CAAst_39d CAAst_39d*
rename CAAst_3d M3410_CAAst_3d CAAst_3d*
rename CAAst_40d M3447_CAAst_40d CAAst_40d*
rename CAAst_41d M3448_CAAst_41d CAAst_41d*
rename CAAst_42d M3449_CAAst_42d CAAst_42d*
rename CAAst_43d M3450_CAAst_43d CAAst_43d*
rename CAAst_44d M3451_CAAst_44d CAAst_44d*
rename CAAst_45d M3452_CAAst_45d CAAst_45d*
rename CAAst_46d M3453_CAAst_46d CAAst_46d*
rename CAAst_47d M3454_CAAst_47d CAAst_47d*
rename CAAst_48d M3455_CAAst_48d CAAst_48d*
rename CAAst_49d M3456_CAAst_49d CAAst_49d*
rename CAAst_4d M3411_CAAst_4d CAAst_4d*
rename CAAst_5d M3412_CAAst_5d CAAst_5d*
rename CAAst_6d M3413_CAAst_6d CAAst_6d*
rename CAAst_7d M3414_CAAst_7d CAAst_7d*
rename CAAst_8d M3415_CAAst_8d CAAst_8d*
rename CAAst_9d M3416_CAAst_9d CAAst_9d*
rename CACst_10d M3468_CACst_10d CACst_10d*
rename CACst_11d M3469_CACst_11d CACst_11d*
rename CACst_12d M3470_CACst_12d CACst_12d*
rename CACst_13d M3471_CACst_13d CACst_13d*
rename CACst_14d M3472_CACst_14d CACst_14d*
rename CACst_15d M3473_CACst_15d CACst_15d*
rename CACst_16d M3474_CACst_16d CACst_16d*
rename CACst_17d M3475_CACst_17d CACst_17d*
rename CACst_18d M3476_CACst_18d CACst_18d*
rename CACst_19d M3477_CACst_19d CACst_19d*
rename CACst_1d M3459_CACst_1d CACst_1d*
rename CACst_20d M3478_CACst_20d CACst_20d*
rename CACst_21d M3479_CACst_21d CACst_21d*
rename CACst_22d M3480_CACst_22d CACst_22d*
rename CACst_23d M3481_CACst_23d CACst_23d*
rename CACst_24d M3482_CACst_24d CACst_24d*
rename CACst_25d M3483_CACst_25d CACst_25d*
rename CACst_26d M3484_CACst_26d CACst_26d*
rename CACst_27d M3485_CACst_27d CACst_27d*
rename CACst_28d M3486_CACst_28d CACst_28d*
rename CACst_29d M3487_CACst_29d CACst_29d*
rename CACst_2d M3460_CACst_2d CACst_2d*
rename CACst_30d M3488_CACst_30d CACst_30d*
rename CACst_31d M3489_CACst_31d CACst_31d*
rename CACst_32d M3490_CACst_32d CACst_32d*
rename CACst_33d M3491_CACst_33d CACst_33d*
rename CACst_34d M3492_CACst_34d CACst_34d*
rename CACst_35d M3493_CACst_35d CACst_35d*
rename CACst_36d M3494_CACst_36d CACst_36d*
rename CACst_37d M3495_CACst_37d CACst_37d*
rename CACst_38d M3496_CACst_38d CACst_38d*
rename CACst_3d M3461_CACst_3d CACst_3d*
rename CACst_4d M3462_CACst_4d CACst_4d*
rename CACst_5d M3463_CACst_5d CACst_5d*
rename CACst_6d M3464_CACst_6d CACst_6d*
rename CACst_7d M3465_CACst_7d CACst_7d*
rename CACst_8d M3466_CACst_8d CACst_8d*
rename CACst_9d M3467_CACst_9d CACst_9d*
rename PL10_B1d M3457_PL10_B1d PL10_B1d*
rename PL10_B2d M3458_PL10_B2d PL10_B2d*
rename CAKst_10d M1401_CAKst_10d CAKst_10d*
rename CAKst_11d M1402_CAKst_11d CAKst_11d*
rename CAKst_12d M1403_CAKst_12d CAKst_12d*
rename CAKst_13d M1404_CAKst_13d CAKst_13d*
rename CAKst_14d M1405_CAKst_14d CAKst_14d*
rename CAKst_15d M1406_CAKst_15d CAKst_15d*
rename CAKst_16d M1407_CAKst_16d CAKst_16d*
rename CAKst_17d M1408_CAKst_17d CAKst_17d*
rename CAKst_18d M1409_CAKst_18d CAKst_18d*
rename CAKst_19d M1410_CAKst_19d CAKst_19d*
rename CAKst_1d M1392_CAKst_1d CAKst_1d*
rename CAKst_20d M1411_CAKst_20d CAKst_20d*
rename CAKst_21d M1412_CAKst_21d CAKst_21d*
rename CAKst_22d M1413_CAKst_22d CAKst_22d*
rename CAKst_23d M1414_CAKst_23d CAKst_23d*
rename CAKst_24d M1415_CAKst_24d CAKst_24d*
rename CAKst_25d M1416_CAKst_25d CAKst_25d*
rename CAKst_26d M1417_CAKst_26d CAKst_26d*
rename CAKst_27d M1418_CAKst_27d CAKst_27d*
rename CAKst_28d M1419_CAKst_28d CAKst_28d*
rename CAKst_29d M1420_CAKst_29d CAKst_29d*
rename CAKst_2d M1393_CAKst_2d CAKst_2d*
rename CAKst_30d M1421_CAKst_30d CAKst_30d*
rename CAKst_31d M1422_CAKst_31d CAKst_31d*
rename CAKst_32d M1423_CAKst_32d CAKst_32d*
rename CAKst_33d M1424_CAKst_33d CAKst_33d*
rename CAKst_34d M1425_CAKst_34d CAKst_34d*
rename CAKst_35d M1426_CAKst_35d CAKst_35d*
rename CAKst_36d M1427_CAKst_36d CAKst_36d*
rename CAKst_37d M1428_CAKst_37d CAKst_37d*
rename CAKst_38d M1429_CAKst_38d CAKst_38d*
rename CAKst_39d M1430_CAKst_39d CAKst_39d*
rename CAKst_3d M1394_CAKst_3d CAKst_3d*
rename CAKst_40d M1431_CAKst_40d CAKst_40d*
rename CAKst_41d M1432_CAKst_41d CAKst_41d*
rename CAKst_42d M1433_CAKst_42d CAKst_42d*
rename CAKst_43d M1434_CAKst_43d CAKst_43d*
rename CAKst_44d M1435_CAKst_44d CAKst_44d*
rename CAKst_45d M1436_CAKst_45d CAKst_45d*
rename CAKst_46d M1437_CAKst_46d CAKst_46d*
rename CAKst_47d M1438_CAKst_47d CAKst_47d*
rename CAKst_48d M1439_CAKst_48d CAKst_48d*
rename CAKst_49d M1440_CAKst_49d CAKst_49d*
rename CAKst_4d M1395_CAKst_4d CAKst_4d*
rename CAKst_50d M1441_CAKst_50d CAKst_50d*
rename CAKst_51d M1442_CAKst_51d CAKst_51d*
rename CAKst_52d M1443_CAKst_52d CAKst_52d*
rename CAKst_5d M1396_CAKst_5d CAKst_5d*
rename CAKst_6d M1397_CAKst_6d CAKst_6d*
rename CAKst_7d M1398_CAKst_7d CAKst_7d*
rename CAKst_8d M1399_CAKst_8d CAKst_8d*
rename CAKst_9d M1400_CAKst_9d CAKst_9d*
rename CAMst_10d M1454_CAMst_10d CAMst_10d*
rename CAMst_11d M1455_CAMst_11d CAMst_11d*
rename CAMst_12d M1456_CAMst_12d CAMst_12d*
rename CAMst_13d M1457_CAMst_13d CAMst_13d*
rename CAMst_14d M1458_CAMst_14d CAMst_14d*
rename CAMst_15d M1459_CAMst_15d CAMst_15d*
rename CAMst_16d M1460_CAMst_16d CAMst_16d*
rename CAMst_17d M1461_CAMst_17d CAMst_17d*
rename CAMst_18d M1462_CAMst_18d CAMst_18d*
rename CAMst_19d M1463_CAMst_19d CAMst_19d*
rename CAMst_1d M1445_CAMst_1d CAMst_1d*
rename CAMst_20d M1464_CAMst_20d CAMst_20d*
rename CAMst_21d M1465_CAMst_21d CAMst_21d*
rename CAMst_22d M1466_CAMst_22d CAMst_22d*
rename CAMst_26d M1467_CAMst_26d CAMst_26d*
rename CAMst_27d M1468_CAMst_27d CAMst_27d*
rename CAMst_28d M1469_CAMst_28d CAMst_28d*
rename CAMst_29d M1470_CAMst_29d CAMst_29d*
rename CAMst_2d M1446_CAMst_2d CAMst_2d*
rename CAMst_30d M1471_CAMst_30d CAMst_30d*
rename CAMst_31d M1472_CAMst_31d CAMst_31d*
rename CAMst_32d M1473_CAMst_32d CAMst_32d*
rename CAMst_33d M1474_CAMst_33d CAMst_33d*
rename CAMst_34d M1475_CAMst_34d CAMst_34d*
rename CAMst_35d M1476_CAMst_35d CAMst_35d*
rename CAMst_36d M1477_CAMst_36d CAMst_36d*
rename CAMst_37d M1478_CAMst_37d CAMst_37d*
rename CAMst_38d M1479_CAMst_38d CAMst_38d*
rename CAMst_39d M1480_CAMst_39d CAMst_39d*
rename CAMst_3d M1447_CAMst_3d CAMst_3d*
rename CAMst_40d M1481_CAMst_40d CAMst_40d*
rename CAMst_41d M1482_CAMst_41d CAMst_41d*
rename CAMst_42d M1483_CAMst_42d CAMst_42d*
rename CAMst_43d M1484_CAMst_43d CAMst_43d*
rename CAMst_44d M1485_CAMst_44d CAMst_44d*
rename CAMst_45d M1486_CAMst_45d CAMst_45d*
rename CAMst_4d M1448_CAMst_4d CAMst_4d*
rename CAMst_5d M1449_CAMst_5d CAMst_5d*
rename CAMst_6d M1450_CAMst_6d CAMst_6d*
rename CAMst_7d M1451_CAMst_7d CAMst_7d*
rename CAMst_8d M1452_CAMst_8d CAMst_8d*
rename CAMst_9d M1453_CAMst_9d CAMst_9d*
rename PL11_B1d M1444_PL11_B1d PL11_B1d*
rename PL11_B2d M1487_PL11_B2d PL11_B2d*
rename CACst_39d M3504_CACst_39d CACst_39d*
rename CACst_40d M3505_CACst_40d CACst_40d*
rename CACst_41d M3506_CACst_41d CACst_41d*
rename CACst_42d M3507_CACst_42d CACst_42d*
rename CACst_43d M3508_CACst_43d CACst_43d*
rename CACst_44d M3509_CACst_44d CACst_44d*
rename CACst_45d M3510_CACst_45d CACst_45d*
rename CACst_46d M3511_CACst_46d CACst_46d*
rename CACst_47d M3512_CACst_47d CACst_47d*
rename CACst_48d M3513_CACst_48d CACst_48d*
rename CACst_49d M3514_CACst_49d CACst_49d*
rename CACst_50d M3515_CACst_50d CACst_50d*
rename CADst_10d M3583_CADst_10d CADst_10d*
rename CADst_11d M3584_CADst_11d CADst_11d*
rename CADst_12d M3585_CADst_12d CADst_12d*
rename CADst_13d M3586_CADst_13d CADst_13d*
rename CADst_14d M3587_CADst_14d CADst_14d*
rename CADst_15d M3588_CADst_15d CADst_15d*
rename CADst_16d M3589_CADst_16d CADst_16d*
rename CADst_17d M3590_CADst_17d CADst_17d*
rename CADst_18d M3591_CADst_18d CADst_18d*
rename CADst_19d M3592_CADst_19d CADst_19d*
rename CADst_1d M3574_CADst_1d CADst_1d*
rename CADst_20d M3593_CADst_20d CADst_20d*
rename CADst_21d M3594_CADst_21d CADst_21d*
rename CADst_22d M3595_CADst_22d CADst_22d*
rename CADst_23d M3596_CADst_23d CADst_23d*
rename CADst_24d M3597_CADst_24d CADst_24d*
rename CADst_25d M3598_CADst_25d CADst_25d*
rename CADst_2d M3575_CADst_2d CADst_2d*
rename CADst_3d M3576_CADst_3d CADst_3d*
rename CADst_4d M3577_CADst_4d CADst_4d*
rename CADst_5d M3578_CADst_5d CADst_5d*
rename CADst_6d M3579_CADst_6d CADst_6d*
rename CADst_7d M3580_CADst_7d CADst_7d*
rename CADst_8d M3581_CADst_8d CADst_8d*
rename CADst_9d M3582_CADst_9d CADst_9d*
rename CAMst_46d M3516_CAMst_46d CAMst_46d*
rename CAMst_47d M3517_CAMst_47d CAMst_47d*
rename CAMst_48d M3518_CAMst_48d CAMst_48d*
rename CAMst_49d M3519_CAMst_49d CAMst_49d*
rename CAMst_50d M3520_CAMst_50d CAMst_50d*
rename CAMst_51d M3521_CAMst_51d CAMst_51d*
rename CANst_10d M3532_CANst_10d CANst_10d*
rename CANst_11d M3533_CANst_11d CANst_11d*
rename CANst_12d M3534_CANst_12d CANst_12d*
rename CANst_13d M3535_CANst_13d CANst_13d*
rename CANst_14d M3536_CANst_14d CANst_14d*
rename CANst_15d M3537_CANst_15d CANst_15d*
rename CANst_16d M3538_CANst_16d CANst_16d*
rename CANst_17d M3539_CANst_17d CANst_17d*
rename CANst_18d M3540_CANst_18d CANst_18d*
rename CANst_19d M3541_CANst_19d CANst_19d*
rename CANst_1d M3523_CANst_1d CANst_1d*
rename CANst_20d M3542_CANst_20d CANst_20d*
rename CANst_21d M3543_CANst_21d CANst_21d*
rename CANst_22d M3544_CANst_22d CANst_22d*
rename CANst_23d M3545_CANst_23d CANst_23d*
rename CANst_24d M3546_CANst_24d CANst_24d*
rename CANst_25d M3547_CANst_25d CANst_25d*
rename CANst_26d M3548_CANst_26d CANst_26d*
rename CANst_27d M3549_CANst_27d CANst_27d*
rename CANst_28d M3550_CANst_28d CANst_28d*
rename CANst_29d M3551_CANst_29d CANst_29d*
rename CANst_2d M3524_CANst_2d CANst_2d*
rename CANst_30d M3552_CANst_30d CANst_30d*
rename CANst_31d M3553_CANst_31d CANst_31d*
rename CANst_32d M3554_CANst_32d CANst_32d*
rename CANst_33d M3555_CANst_33d CANst_33d*
rename CANst_34d M3556_CANst_34d CANst_34d*
rename CANst_35d M3557_CANst_35d CANst_35d*
rename CANst_36d M3558_CANst_36d CANst_36d*
rename CANst_37d M3559_CANst_37d CANst_37d*
rename CANst_38d M3560_CANst_38d CANst_38d*
rename CANst_39d M3561_CANst_39d CANst_39d*
rename CANst_3d M3525_CANst_3d CANst_3d*
rename CANst_40d M3562_CANst_40d CANst_40d*
rename CANst_41d M3563_CANst_41d CANst_41d*
rename CANst_42d M3564_CANst_42d CANst_42d*
rename CANst_43d M3565_CANst_43d CANst_43d*
rename CANst_44d M3566_CANst_44d CANst_44d*
rename CANst_45d M3567_CANst_45d CANst_45d*
rename CANst_46d M3568_CANst_46d CANst_46d*
rename CANst_47d M3569_CANst_47d CANst_47d*
rename CANst_48d M3570_CANst_48d CANst_48d*
rename CANst_49d M3571_CANst_49d CANst_49d*
rename CANst_4d M3526_CANst_4d CANst_4d*
rename CANst_50d M3572_CANst_50d CANst_50d*
rename CANst_5d M3527_CANst_5d CANst_5d*
rename CANst_6d M3529_CANst_6d CANst_6d*
rename CANst_7d M3528_CANst_7d CANst_7d*
rename CANst_8d M3530_CANst_8d CANst_8d*
rename CANst_9d M3531_CANst_9d CANst_9d*
rename PL12_B1d M3522_PL12_B1d PL12_B1d*
rename PL12_B2d M3573_PL12_B2d PL12_B2d*
rename PL12_B3d M3599_PL12_B3d PL12_B3d*
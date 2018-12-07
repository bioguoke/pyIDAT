#coding=utf-8
#!/usr/bin/python
# Created on 
# Author: Zhihua Pei
# Organization: www.bioguoke.com
# website: www.zilhua.com
# github: github.com/zilhua

# src : https://github.com/HenrikBengtsson/illuminaio
# R to python


import struct
import os,sys
import pandas as pd

class pyIDAT(object):
    fileSize, versionNumber, nFields, fields, nSNPsRead, Quants, \
        MidBlock, RedGreen, Barcode, ChipType, RunInfo, Unknowns = [None] * 12

def readIDAT(file=None, what="all",**kwargs):
    # what = "all", "IlluminaID", "nSNPsRead"
    # file = u"xxxxxxx.idat"
    with open(file, "rb") as f:
        version = readBin(f, n=1)
    if version[0] == 3:
        # 非加密IDAT文件
        return readIDAT_nonenc(file, what=what)
    else:
        # 其他格式的文件
        raise IOError("err: no idat fotmat!")

def readIDAT_nonenc(file, what="all"):
    fields = []
    colname = ["fieldCode", "byteOffset", "Bytes"]
    rowname = []
    knownCodes = {
        1000  : "nSNPsRead",
        102   : "IlluminaID",
        103   : "SD",
        104   : "Mean",
        107   : "NBeads",
        200   : "MidBlock",
        300   : "RunInfo",
        400   : "RedGreen",
        401   : "MostlyNull",# 'Manifest', cf [1].
        402   : "Barcode",
        403   : "ChipType",
        404   : "MostlyA",   # 'Stripe', cf [1].
        405   : "Unknown.1",
        406   : "Unknown.2", # 'Sample ID', cf [1].
        407   : "Unknown.3",
        408   : "Unknown.4", # 'Plate', cf [1].
        409   : "Unknown.5", # 'Well', cf [1].
        410   : "Unknown.6",
        510   : "Unknown.7",
    }
    idat_data = pyIDAT()
    idat_data.fileSize = os.path.getsize(file)
    with open(file, "rb") as f:
        version = readBin(f, n=1)
        version = version[0]
        if version != 3:
            raise IOError("Err: Cant read idat file")
        # 获取字段长度
        nFields = _readInt(f, n=1)[0]
        nNewFields = 1
        for i in xrange(nFields):
            fieldCode  = _readShort(f, n=1)[0]
            byteOffset = _readLong(f, n=1)[0]
            fields.append([fieldCode, byteOffset, 0])
            if fieldCode in knownCodes:
                rowname.append(knownCodes[fieldCode])
            else:
                rowname.append("newField.{0}".format(nNewFields))
                nNewFields += 1
        fields = pd.DataFrame(data=fields, columns=colname, index=rowname)
        if fields["byteOffset"]["nSNPsRead"] != min(fields["byteOffset"]):
            raise IOError("read idat file err")
        f.seek(fields["byteOffset"]["nSNPsRead"], 0)
        nSNPsRead = _readInt(f, n=1)[0]
        if what == "nSNPsRead":
            return nSNPsRead
        if what == "IlluminaID":
            f.seek(fields["byteOffset"]["IlluminaID"], 0)
            res = _readField(f, "IlluminaID", nSNPsRead)
            return res
        def _find_field(fieldname):
            pos = fields["byteOffset"][fieldname]
            f.seek(pos, 0)
            return _readField(f, fieldname, nSNPsRead)
        rowname_sorted = fields.sort_values(by="byteOffset").index.tolist()
        res = { k: _find_field(k) for k in rowname_sorted if k != "nSNPsRead"}
        Unknowns = {
            "MostlyNull"  : res["MostlyNull"],
            "MostlyA"     : res["MostlyA"],
            "Unknown.1"   : res["Unknown.1"],
            "Unknown.2"   : res["Unknown.2"],
            "Unknown.3"   : res["Unknown.3"],
            "Unknown.4"   : res["Unknown.4"],
            "Unknown.5"   : res["Unknown.5"],
        }
        Quants = pd.DataFrame(
            {"Mean" : res["Mean"], "SD" : res["SD"], "NBeads" : res["NBeads"]},
            index=map(str,res["IlluminaID"]))
        idat_data.versionNumber = version
        idat_data.nFields       = nFields
        idat_data.fields        = fields
        idat_data.nSNPsRead     = nSNPsRead
        idat_data.Quants        = Quants
        idat_data.MidBlock      = res["MidBlock"]
        idat_data.RedGreen      = res["RedGreen"]
        idat_data.Barcode       = res["Barcode"]
        idat_data.ChipType      = res["ChipType"]
        idat_data.RunInfo       = res["RunInfo"]
        idat_data.Unknowns      = Unknowns
    return idat_data

def readBin(f, size=4, n=1, signed ="i", endian = "little"):
    r = []
    count = 0
    byte = f.read(4)
    if byte != "IDAT":
        raise IOError("Cannot read IDAT file. File format error. Unknown magic: {0}".format(byte))
    while count < n and byte != b"":
        byte = f.read(size)
        if not byte:
            return r + [0] * (n - len(r))
        i = struct.unpack("<i", byte)[0]
        r.append(i)
        count += 1
        byte = f.read(size)
    return r

def _readLong(con, n=1):
    r, count = [], 0
    byte = True
    while count < n and byte != b"":
        byte = con.read(4)
        i = struct.unpack("<i", byte)[0]
        r.append(i)
        count += 1
        byte = con.read(4)
    return r

def _readInt(con, n=1):
    r, count = [], 0
    byte = True
    while count < n and byte != b"":
        byte = con.read(4)
        i = struct.unpack("<i", byte)[0]
        r.append(i)
        count += 1
    return r

def _readShort(con, n=1):
    r, count = [], 0
    byte = True
    while count < n and byte != b"":
        byte = con.read(2)
        i = struct.unpack("<H", byte)[0]
        r.append(i)
        count += 1
    return r

def _readByte(con, n=1):
    r, count = [], 0
    byte = True
    while count < n and byte != b"":
        byte = con.read(1)
        if not byte:
            return r + [0] * (n-len(r))
        i = struct.unpack("<B", byte)[0]
        r.append(i)
        count += 1
        # byte = con.read(1)
    return r

def _readString(con):
    m = _readByte(con, n=1)[0]
    n = m % 128
    shift = 0
    while(m / 128 == 1 ):
        m = _readByte(con, n=1)[0]
        shift = shift + 7
        k = (m % 128) * (2 ** shift)
        n = n + k
    res = con.read(n)
    return res

def _readField(con, field, nSNPsRead=1):
    if field == "RunInfo":
        nRunInfoBlocks = _readInt(con, n=1)[0]
        RunInfo = []
        for i in xrange(nRunInfoBlocks):
            RunInfo.append(
                [_readString(con) for j in xrange(5)]
            )
        RunInfo = pd.DataFrame(data=RunInfo, columns=["RunTime", "BlockType", "BlockPars", "BlockCode", "CodeVersion"])
        return RunInfo
    dic = {
        "IlluminaID" : lambda x:_readInt(x, n=nSNPsRead),
        "SD"         : lambda x:_readShort(x, n=nSNPsRead),
        "Mean"       : lambda x:_readShort(x, n=nSNPsRead),
        "NBeads"     : lambda x:_readByte(x, n=nSNPsRead),
        "MidBlock"   : lambda x:_readInt(x, n=_readInt(x, n=1)[0]),
        "RedGreen"   : lambda x:_readInt(x, n=1)[0],
        "MostlyNull" : lambda x:_readString(x),
        "Barcode"    : lambda x:_readString(x),
        "ChipType"   : lambda x:_readString(x),
        "MostlyA"    : lambda x:_readString(x),
        "Unknown.1"  : lambda x:_readString(x),
        "Unknown.2"  : lambda x:_readString(x),
        "Unknown.3"  : lambda x:_readString(x),
        "Unknown.4"  : lambda x:_readString(x),
        "Unknown.5"  : lambda x:_readString(x),
        "Unknown.6"  : lambda x:_readString(x),
        "Unknown.7"  : lambda x:_readString(x),
    }
    try:
        return dic.get(field)(con)
    except Exception,e:
        raise IOError("no field:".format(field))




if __name__ == "__main__":
    idat_data = readIDAT(what="all")
    # fileSize, versionNumber, nFields, fields, nSNPsRead, Quants, \
    #         MidBlock, RedGreen, Barcode, ChipType, RunInfo, Unknowns = [None] * 12
    print idat_data.RunInfo
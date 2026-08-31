# -*- coding: utf-8 -*-
"""生成 TerraPos 应用图标(固定路径, 无命令行参数)。

设计: 深青底圆角方, 山体按输出图例配色分三带(坡上橙/坡中红棕/坡下深棕),
山顶白色峰顶种子点, 底部盆地蓝水带 —— 直接呼应产品功能(8 类地形部位)。
输出: rust/topo_app/assets/icon.ico + icon-256.png
"""
import os
from PIL import Image, ImageDraw

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "rust", "topo_app", "assets")

S = 256
# 产品输出图例配色(terrain_position cmap)
C_TOP = (250, 165, 60)     # 山地坡上 橙
C_MID = (222, 100, 50)     # 山地坡中 红棕
C_BOT = (150, 118, 89)     # 山地坡下 棕(略提亮避免过暗)
C_BASIN = (51, 178, 229)   # 山间盆地 蓝
C_BG_TOP = (31, 59, 87)    # 背景深青上
C_BG_BOT = (24, 42, 62)    # 背景深青下
C_PEAK = (255, 255, 255)


def main():
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)

    # 背景: 圆角方 + 竖直渐变
    grad = Image.new("RGBA", (S, S))
    gd = ImageDraw.Draw(grad)
    for y in range(S):
        t = y / (S - 1)
        c = tuple(int(C_BG_TOP[i] + (C_BG_BOT[i] - C_BG_TOP[i]) * t) for i in range(3)) + (255,)
        gd.line([(0, y), (S, y)], fill=c)
    mask = Image.new("L", (S, S), 0)
    ImageDraw.Draw(mask).rounded_rectangle([4, 4, S - 4, S - 4], radius=52, fill=255)
    img.paste(grad, (0, 0), mask)

    # 山体三带(等高线式水平分带): 三角峰 (128,52) → 底 (34,212)-(222,212)
    peak_y, base_y = 52, 212
    xl, xr = 34, 222

    def edges(y):
        w = (y - peak_y) / (base_y - peak_y)
        return (128 + (xl - 128) * w, 128 + (xr - 128) * w)

    def band(y0, y1, color):
        a0, b0 = edges(y0)
        a1, b1 = edges(y1)
        d.polygon([(a0, y0), (b0, y0), (b1, y1), (a1, y1)], fill=color)

    band(peak_y, 104, C_TOP)
    band(104, 158, C_MID)
    band(158, base_y, C_BOT)

    def mix(c, alpha):
        """白色按 alpha 比例预混合到实色(规避 ImageDraw 不做 alpha 混合)。"""
        return tuple(int(c[i] + (255 - c[i]) * alpha) for i in range(3))

    # 带间等高线分隔(白色预混合实色)
    for y in (104, 158):
        a, b = edges(y)
        d.line([(a, y), (b, y)], fill=(215, 218, 222), width=3)

    # 左坡受光面提亮(预混合, 不用半透明层)
    hl = ImageDraw.Draw(img)
    for yy in range(peak_y + 2, base_y, 2):
        a, b = edges(yy)
        xa = a + (b - a) * 0.18
        xb = a + (b - a) * 0.42
        base = C_TOP if yy < 104 else (C_MID if yy < 158 else C_BOT)
        hl.line([(xa, yy), (xb, yy)], fill=mix(base, 0.16), width=2)

    # 峰顶种子点(白心橙环)
    d.ellipse([118, 42, 138, 62], fill=C_PEAK, outline=(240, 120, 40, 255), width=4)

    # 底部盆地蓝水带 + 两条波纹
    d.rounded_rectangle([24, 222, 232, 236], radius=7, fill=C_BASIN + (235,))
    for i, y in enumerate((227, 231)):
        a = 36 + i * 10
        b = 220 - i * 10
        d.line([(a, y), (b, y)], fill=(255, 255, 255, 110), width=2)

    os.makedirs(OUT_DIR, exist_ok=True)
    png = os.path.join(OUT_DIR, "icon-256.png")
    ico = os.path.join(OUT_DIR, "icon.ico")
    img.save(png)
    img.save(ico, sizes=[(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (16, 16)])
    print("已生成:", os.path.abspath(png))
    print("已生成:", os.path.abspath(ico))


if __name__ == "__main__":
    main()

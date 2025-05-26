import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class CatenaryLumpedMass:
    """三本係留索のランプドマス法による初期位置計算クラス"""
    
    def __init__(self, num_segments=20):
        self.num_segments = num_segments  # セグメント数（内部ノード数 = num_segments - 1）
        self.gravity = 9.80665  # 重力加速度 [m/s²]
        
    def myasinh(self, val):
        """逆双曲線正弦関数"""
        return math.log(val + math.sqrt(val * val + 1.0))
    
    def myacosh(self, val):
        """逆双曲線余弦関数"""
        if val < 1.0:
            val = 1.0  # 引数保護
        return math.log(val + math.sqrt(val + 1.0) * math.sqrt(val - 1.0))
    
    def myatanh(self, val):
        """逆双曲線正接関数"""
        return 0.5 * math.log((1.0 + val) / (1.0 - val))
    
    def funcd_catenary(self, x_param, xacc, d_geom, l_geom):
        """カテナリー方程式とその導関数を計算"""
        max_internal_iter = 1000
        p0_angle = 0.0
        
        # 水平張力が無い状態
        if x_param == 0.0:
            f_val = -d_geom
            df_val = 0.0
            return f_val, df_val, p0_angle
        
        # 水平張力がある状態
        elif x_param > 0.0:
            # 特殊なケース：l_geom <= 0の場合
            if l_geom <= 0.0:
                X_1_internal = 1.0/x_param + 1.0
                if X_1_internal < 1.0:
                    X_1_internal = 1.0
                
                term_sqrt1 = 1.0 + 2.0*x_param
                if term_sqrt1 < 0.0:
                    term_sqrt1 = 0.0
                
                term_sqrt2 = X_1_internal*X_1_internal - 1.0
                if term_sqrt2 < 0.0:
                    term_sqrt2 = 0.0
                
                f_val = x_param*self.myacosh(X_1_internal) - math.sqrt(term_sqrt1) + 1.0 - d_geom
                
                if abs(x_param * math.sqrt(term_sqrt2)) < 1e-12 or abs(math.sqrt(term_sqrt1)) < 1e-12:
                    df_val = self.myacosh(X_1_internal)
                else:
                    df_val = self.myacosh(X_1_internal) - 1.0/math.sqrt(term_sqrt1) - 1.0/(x_param*math.sqrt(term_sqrt2))
                
                p0_angle = 0.0
            
            # 一般的なケース
            else:
                # 海底との接触がある可能性がある場合
                if x_param > (l_geom*l_geom - 1.0) / 2.0:
                    p0_angle = 0.0
                    for i in range(max_internal_iter):
                        cos_p0 = math.cos(p0_angle)
                        if abs(cos_p0) < 1e-9:
                            df1_internal = 1.0
                            break
                        
                        func1_internal = 1.0/x_param + 1.0/cos_p0
                        term_in_sqrt_f1 = func1_internal*func1_internal - 1.0
                        if term_in_sqrt_f1 < 0.0:
                            term_in_sqrt_f1 = 0.0
                        
                        f1_internal = x_param*(math.sqrt(term_in_sqrt_f1) - math.tan(p0_angle)) - l_geom
                        
                        if abs(cos_p0) < 1e-9 or term_in_sqrt_f1 < 1e-12:
                            df1_internal = 1.0
                            break
                        
                        df1_internal = x_param * (func1_internal * math.tan(p0_angle) / (cos_p0 * math.sqrt(term_in_sqrt_f1)) - (math.tan(p0_angle)*math.tan(p0_angle)) - 1.0)
                        
                        if abs(df1_internal) < 1e-9:
                            break
                        
                        p0_angle = p0_angle - f1_internal/df1_internal
                        
                        cos_p0 = math.cos(p0_angle)
                        if abs(cos_p0) < 1e-9:
                            break
                        
                        func1_internal = 1.0/x_param + 1.0/cos_p0
                        term_in_sqrt_f1 = func1_internal*func1_internal - 1.0
                        if term_in_sqrt_f1 < 0.0:
                            term_in_sqrt_f1 = 0.0
                        
                        f1_internal = x_param*(math.sqrt(term_in_sqrt_f1) - math.tan(p0_angle)) - l_geom
                        
                        if abs(f1_internal) < xacc:
                            break
                    
                    X_2_internal = l_geom/x_param + math.tan(p0_angle)
                    X_3_internal = math.tan(p0_angle)
                    f_val = x_param*(self.myasinh(X_2_internal) - self.myasinh(X_3_internal)) - l_geom + 1.0 - d_geom
                    
                    term_in_sqrt_df = X_2_internal*X_2_internal + 1.0
                    if abs(x_param * math.sqrt(term_in_sqrt_df)) < 1e-12:
                        df_val = 1.0
                    else:
                        df_val = self.myasinh(X_2_internal) - self.myasinh(X_3_internal) - l_geom/(x_param*math.sqrt(term_in_sqrt_df))
                
                # 海底との接触がない場合：単純なカテナリー
                else:
                    X_5_internal = 1.0/x_param + 1.0
                    if X_5_internal < 1.0:
                        X_5_internal = 1.0
                    
                    term_sqrt1 = 1.0 + 2.0*x_param
                    if term_sqrt1 < 0.0:
                        term_sqrt1 = 0.0
                    
                    term_sqrt2 = X_5_internal*X_5_internal - 1.0
                    if term_sqrt2 < 0.0:
                        term_sqrt2 = 0.0
                    
                    f_val = x_param*self.myacosh(X_5_internal) - math.sqrt(term_sqrt1) + 1.0 - d_geom
                    
                    if abs(x_param * math.sqrt(term_sqrt2)) < 1e-12 or abs(math.sqrt(term_sqrt1)) < 1e-12:
                        df_val = self.myacosh(X_5_internal)
                    else:
                        df_val = self.myacosh(X_5_internal) - 1.0/math.sqrt(term_sqrt1) - 1.0/(x_param*math.sqrt(term_sqrt2))
                    
                    p0_angle = 0.0
        
        else:
            raise ValueError("x_param must be non-negative")
        
        return f_val, df_val, p0_angle
    
    def rtsafe_catenary(self, x1_bounds, x2_bounds, xacc_tol, d_geom, l_geom):
        """Newton-Raphson法による根の探索"""
        MAXIT = 1000
        
        # 境界での関数値
        fl, _, p1 = self.funcd_catenary(x1_bounds, xacc_tol, d_geom, l_geom)
        fh, _, p2 = self.funcd_catenary(x2_bounds, xacc_tol, d_geom, l_geom)
        
        # 同符号ならエラー
        if (fl > 0.0 and fh > 0.0) or (fl < 0.0 and fh < 0.0):
            raise ValueError("Root must be bracketed in rtsafe")
        
        # どちらかが0なら答え
        if fl == 0.0:
            return x1_bounds, p1
        if fh == 0.0:
            return x2_bounds, p2
        
        # fl < 0となるようにxlとxhを設定
        if fl < 0.0:
            xl = x1_bounds
            xh = x2_bounds
        else:
            xh = x1_bounds
            xl = x2_bounds
        
        # 中点を初期推定値とする
        rts = 0.5*(x1_bounds + x2_bounds)
        dxold = abs(x2_bounds - x1_bounds)
        dx = dxold
        f, df, p0_angle_out = self.funcd_catenary(rts, xacc_tol, d_geom, l_geom)
        
        for j in range(MAXIT):
            if (((rts - xh)*df - f)*((rts - xl)*df - f) > 0.0) or (abs(2.0*f) > abs(dxold*df)):
                dxold = dx
                dx = 0.5*(xh - xl)
                rts = xl + dx
                if xl == rts:
                    return rts, p0_angle_out
            else:
                dxold = dx
                if abs(df) < 1e-12:
                    dx = 0.5*(xh - xl)
                    rts = xl + dx
                    if xl == rts:
                        return rts, p0_angle_out
                else:
                    dx = f/df
                
                temp = rts
                rts -= dx
                if temp == rts:
                    return rts, p0_angle_out
            
            if abs(dx) < xacc_tol:
                return rts, p0_angle_out
            
            f, df, p0_angle_out = self.funcd_catenary(rts, xacc_tol, d_geom, l_geom)
            
            if f < 0.0:
                xl = rts
            else:
                xh = rts
        
        raise ValueError("Maximum iterations exceeded in rtsafe")
    
    def calculate_single_mooring_positions(self, fairlead_pos, anchor_pos, total_length, unit_weight, rtsafe_accuracy):
        """単一係留索の初期位置を計算"""
        
        # アンカーを原点とした座標系でのフェアリーダー位置
        FP_AP = np.array(fairlead_pos) - np.array(anchor_pos)
        
        # 鉛直距離と水平距離
        h = abs(FP_AP[2])  # Z方向の距離
        L_APFP = math.sqrt(FP_AP[0]**2 + FP_AP[1]**2)  # 水平距離
        
        # 水平面での単位ベクトル
        if L_APFP > 1e-12:
            horizontal_dir = np.array([FP_AP[0] / L_APFP, FP_AP[1] / L_APFP, 0.0])
        else:
            horizontal_dir = np.array([0.0, 0.0, 0.0])
        
        # カテナリー理論パラメータ
        L0_APFP = total_length - h  # 水平張力が0になるときの水平距離
        delta = L_APFP - L0_APFP
        
        if h > 1e-12:
            d = delta / h
            l = total_length / h
        else:
            d = 0.0
            l = 0.0
        
        # 水平張力パラメータの計算
        x_param = 0.0
        p0 = 0.0
        H = 0.0  # 水平張力
        
        if d <= 0.0:
            # たるんでいる場合
            H = 0.0
        elif h > 1e-12 and d < (math.sqrt(l*l - 1.0) - (l - 1.0)):
            # カテナリー形状
            x1 = 0.0
            x2 = 1.0e6
            x_param, p0 = self.rtsafe_catenary(x1, x2, rtsafe_accuracy, d, l)
            H = x_param * unit_weight * h
        
        # 各ノードの位置計算
        node_positions = []
        node_positions.append(fairlead_pos)  # フェアリーダー位置
        
        if H > 1e-12:
            # カテナリー形状に沿って配置
            a = H / unit_weight  # カテナリーパラメータ
            
            # APからの弧長を計算
            arc_length = [total_length * i / self.num_segments for i in range(self.num_segments + 1)]
            
            # 弧長に対応する位置の計算（内部ノード：1からnum_segments-1まで）
            for i in range(1, self.num_segments):
                s = total_length - arc_length[i]  # FPからの弧長
                
                # カテナリー曲線に沿った位置の計算
                if p0 < 1e-12:
                    # 海底接触無しの場合
                    beta = s / a
                    if beta < 50.0:
                        x_local = a * self.myasinh(math.sinh(L_APFP / a) - math.sinh(beta))
                        z_local = a * (math.cosh((L_APFP - x_local) / a) - math.cosh(L_APFP / a))
                    else:
                        x_local = L_APFP - s
                        z_local = 0.0
                else:
                    # 海底接触ありの場合（簡略化）
                    ratio = i / self.num_segments
                    x_local = L_APFP * (1.0 - ratio)
                    z_local = h * (1.0 - ratio)
                
                # グローバル座標系における各ノードの位置
                position = np.array(anchor_pos) + horizontal_dir * x_local + np.array([0.0, 0.0, z_local])
                node_positions.append(position.tolist())
        
        else:
            # 垂直に垂れ下がる場合
            for i in range(1, self.num_segments):
                ratio = i / self.num_segments
                z_offset = -h * ratio
                position = np.array(fairlead_pos) + np.array([0.0, 0.0, z_offset])
                node_positions.append(position.tolist())
        
        return node_positions
    
    def calculate_three_mooring_positions(self, fairlead_positions, anchor_positions, 
                                        total_lengths, unit_weights, rtsafe_accuracy):
        """三本係留索の初期位置を計算"""
        
        if len(fairlead_positions) != 3 or len(anchor_positions) != 3:
            raise ValueError("Three fairlead and anchor positions must be provided")
        
        if len(total_lengths) != 3 or len(unit_weights) != 3:
            raise ValueError("Three total lengths and unit weights must be provided")
        
        all_mooring_positions = []
        
        for i in range(3):
            print(f"\n=== 係留索 {i+1} の計算 ===")
            print(f"フェアリーダー位置: {fairlead_positions[i]}")
            print(f"アンカー位置: {anchor_positions[i]}")
            print(f"全長: {total_lengths[i]} m")
            print(f"単位重量: {unit_weights[i]} N/m")
            
            positions = self.calculate_single_mooring_positions(
                fairlead_positions[i], 
                anchor_positions[i], 
                total_lengths[i], 
                unit_weights[i], 
                rtsafe_accuracy
            )
            
            all_mooring_positions.append(positions)
            print(f"内部ノード数: {len(positions) - 1}")
        
        return all_mooring_positions
    
    def visualize_moorings(self, all_mooring_positions, fairlead_positions, anchor_positions):
        """三本係留索の3D可視化"""
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        colors = ['red', 'blue', 'green']
        
        for i, positions in enumerate(all_mooring_positions):
            # 係留索の線
            x_coords = [pos[0] for pos in positions]
            y_coords = [pos[1] for pos in positions]
            z_coords = [pos[2] for pos in positions]
            
            # アンカーも含めて描画
            x_coords.append(anchor_positions[i][0])
            y_coords.append(anchor_positions[i][1])
            z_coords.append(anchor_positions[i][2])
            
            ax.plot(x_coords, y_coords, z_coords, 
                   color=colors[i], linewidth=2, label=f'係留索 {i+1}')
            
            # 内部ノード
            ax.scatter(x_coords[1:-1], y_coords[1:-1], z_coords[1:-1], 
                      color=colors[i], s=30, alpha=0.6)
            
            # フェアリーダー
            ax.scatter(fairlead_positions[i][0], fairlead_positions[i][1], fairlead_positions[i][2], 
                      color=colors[i], s=100, marker='s', label=f'FP {i+1}')
            
            # アンカー
            ax.scatter(anchor_positions[i][0], anchor_positions[i][1], anchor_positions[i][2], 
                      color=colors[i], s=100, marker='^', label=f'AP {i+1}')
        
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
        ax.set_title('三本係留索の初期配置')
        ax.legend()
        ax.grid(True)
        
        plt.tight_layout()
        plt.show()

def main():
    """メイン関数"""
    
    # パラメータ設定
    catenary = CatenaryLumpedMass(num_segments=20)  # 20セグメント（19個の内部ノード）
    
    # 入力パラメータ
    print("=== 三本係留索ランプドマス法初期位置計算 ===\n")
    
    # 係留索パラメータ
    total_lengths = [902.2, 902.2, 902.2]  # 全長 [m]
    line_mass = 77.71
    rho = 1025
    line_diameter = 0.09017
    gravity = 9.80665
    line_weight = (line_mass - rho*0.25*math.pi*line_diameter^2)*gravity
    unit_weights = [line_weight, line_weight, line_weight]   # 単位重量 [N/m]
    rtsafe_accuracy = 1e-4  # 計算精度

    anchor_radius = 853.87
    anchor_depth = -320
    fairleader_radius = 5.2
    fairleader_depth = -70.0

    xap1 = anchor_radius
    yap1 = 0.0
    zap1 = anchor_depth

    xap2 = anchor_radius*math.cos(math.pi * 2.0/3.0)
    yap2 = anchor_radius*math.sin(math.pi * 2.0/3.0)
    zap2 = anchor_depth

    xap3 = anchor_radius*math.cos(-math.pi * 2.0/3.0)
    yap3 = anchor_radius*math.sin(-math.pi * 2.0/3.0)
    zap3 = anchor_depth

    xfp1 = fairleader_radius
    yfp1 = 0.0
    zfp1 = fairleader_depth

    xfp2 = fairleader_radius*math.cos(math.pi * 2.0/3.0)
    yfp2 = fairleader_radius*math.sin(math.pi * 2.0/3.0)
    zfp2 = fairleader_depth

    xfp3 = fairleader_radius*math.cos(-math.pi * 2.0/3.0)
    yfp3 = fairleader_radius*math.sin(-math.pi * 2.0/3.0)
    zfp3 = fairleader_depth
    
    # フェアリーダー位置（三次元）
    fairlead_positions = [
        [xfp1, yfp1, zfp1],    # FP1
        [xfp2, yfp2, zfp2],    # FP2
        [xfp3, yfp3, zfp3]     # FP3
    ]
    
    # アンカー位置（三次元）
    anchor_positions = [
        [xap1, yap1, zap1],    # AP1
        [xap2, yap2, zap2],    # AP2
        [xap3, yap3, zap3]     # AP3
    ]
    
    print("入力パラメータ:")
    print(f"セグメント数: {catenary.num_segments}")
    print(f"内部ノード数（一本あたり）: {catenary.num_segments - 1}")
    print(f"全長: {total_lengths} m")
    print(f"単位重量: {unit_weights} N/m")
    print(f"計算精度: {rtsafe_accuracy}")
    print(f"フェアリーダー位置: {fairlead_positions}")
    print(f"アンカー位置: {anchor_positions}")
    
    # 三本係留索の初期位置計算
    try:
        all_positions = catenary.calculate_three_mooring_positions(
            fairlead_positions, anchor_positions, 
            total_lengths, unit_weights, rtsafe_accuracy
        )
        
        # 結果出力
        print("\n=== 計算結果 ===")
        for i, positions in enumerate(all_positions):
            print(f"\n係留索 {i+1} の内部ノード位置:")
            print("ノード番号 | X座標 [m] | Y座標 [m] | Z座標 [m]")
            print("-" * 50)
            
            # フェアリーダー
            fp = positions[0]
            print(f"FP      | {fp[0]:8.3f} | {fp[1]:8.3f} | {fp[2]:8.3f}")
            
            # 内部ノード
            for j in range(1, len(positions)):
                pos = positions[j]
                print(f"Node{j:2d} | {pos[0]:8.3f} | {pos[1]:8.3f} | {pos[2]:8.3f}")
            
            # アンカー
            ap = anchor_positions[i]
            print(f"AP      | {ap[0]:8.3f} | {ap[1]:8.3f} | {ap[2]:8.3f}")
        
        # 可視化
        catenary.visualize_moorings(all_positions, fairlead_positions, anchor_positions)
        
        # CSVファイルに保存
        import csv
        with open('mooring_positions.csv', 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['係留索番号', 'ノード種別', 'ノード番号', 'X座標[m]', 'Y座標[m]', 'Z座標[m]'])
            
            for i, positions in enumerate(all_positions):
                # フェアリーダー
                fp = positions[0]
                writer.writerow([i+1, 'Fairlead', 0, fp[0], fp[1], fp[2]])
                
                # 内部ノード
                for j in range(1, len(positions)):
                    pos = positions[j]
                    writer.writerow([i+1, 'Internal', j, pos[0], pos[1], pos[2]])
                
                # アンカー
                ap = anchor_positions[i]
                writer.writerow([i+1, 'Anchor', -1, ap[0], ap[1], ap[2]])
        
        print(f"\n結果を 'mooring_positions.csv' に保存しました。")
        
    except Exception as e:
        print(f"エラーが発生しました: {e}")

if __name__ == "__main__":
    main()

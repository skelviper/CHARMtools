import matplotlib.pyplot as plt
import seaborn as sns
import math
from typing import Optional, Sequence, Tuple, Dict
import numpy as np
import pandas as pd
from scipy import sparse
from pandas.api.types import CategoricalDtype

# 从utils导入通用函数
from ..utils.helper import parse_genome_coord, natural_chr_key

def _get_X(adata, layer):
    """获取 AnnData 的 X 或指定 layer"""
    if layer is None:
        return adata.X
    return adata.layers[layer]

class MultiCell3DVisualization:
    """
    Visualization methods for MultiCell3D class.
    """
    
    def plot_diff(self, diff, chrom_plot):
        """
        Plot the difference between two groups of cells.
        """
        df = diff.query('chrom == @chrom_plot')
        plt.figure(figsize=(6, 4))
        plt.plot(df['pos'], df['mean_group1'], label='Group 1', color='blue', alpha=1)
        plt.plot(df['pos'], df['mean_group2'], label='Group 2', color='red', alpha=1)

        significant_points = df[df['p_value_adj'] < 0.05]
        for index, row in significant_points.iterrows():
            plt.axvspan(row['pos'], row['pos'] + 1000000, color='grey', alpha=0.2)

        plt.ylim(0.7, 1.3)
        plt.xlabel(chrom_plot)
        plt.ylabel('Radial position')
        plt.legend()

        plt.show()
    
    def plot_radial_position(
        self,
        adata_name: str,
        group_by: str,
        facet_by: str = "chrom",
        genome_coord: str = "",
        *,
        layer: Optional[str] = None,
        groups: Optional[Sequence[str]] = None,
        color_palette: Optional[Dict[str, str]] = None,
        groups_as_rows: bool = False,
        merge_alleles: bool = True,
        facets: Optional[Sequence[str]] = None,
        figsize: Tuple[float, float] = (6, 4),
        dpi: int = 150,
        alpha: float = 1,
        s: float = 3.0,
        legend_outside: bool = True,
        ylim: Optional[Tuple[float, float]] = None,
    ):
        """
        绘制径向位置数据的高级可视化函数
        
        Parameters:
        -----------
        adata_name : str
            AnnData矩阵的名称，将从对象的adata_dict中获取
        group_by : str
            用于分组的obs列名
        facet_by : str, default "chrom"
            用于分面的var列名（通常是染色体）
        genome_coord : str, default ""
            基因组坐标范围，格式如"chr1:1000000-2000000"
        layer : str, optional
            使用的数据层，默认使用X
        groups : Sequence[str], optional
            要绘制的特定分组
        color_palette : Dict[str, str], optional
            分组颜色映射
        groups_as_rows : bool, default False
            是否将分组作为行显示
        merge_alleles : bool, default True
            是否合并等位基因
        facets : Sequence[str], optional
            要绘制的特定染色体
        figsize : Tuple[float, float], default (6, 4)
            图形大小
        dpi : int, default 150
            图形分辨率
        alpha : float, default 1
            点的透明度
        s : float, default 3.0
            点的大小
        legend_outside : bool, default True
            是否将图例放在外侧
        ylim : Tuple[float, float], optional
            y轴的显示范围，格式为(下限, 上限)。如果为None，则自动调整
            
        Returns:
        --------
        matplotlib.figure.Figure
            创建的图形对象
        """
        
        # 获取 AnnData 对象
        adata = self.get_anndata(adata_name)
        
        # 1. 数据选择：基因组坐标过滤
        coord_info = parse_genome_coord(genome_coord)
        if coord_info:
            target_chrom, start_pos, end_pos = coord_info
            # 过滤 var：染色体匹配 + 位置范围
            chrom_mask = adata.var[facet_by] == target_chrom
            if "pos" in adata.var.columns:
                pos_mask = (adata.var["pos"] >= start_pos) & (adata.var["pos"] <= end_pos)
                var_mask = chrom_mask & pos_mask
            else:
                var_mask = chrom_mask
            adata_sub = adata[:, var_mask]
        else:
            adata_sub = adata
        
        if adata_sub.n_vars == 0:
            raise ValueError(f"No variables found for genome coordinate: {genome_coord}")
        
        # 2. 等位基因合并
        if merge_alleles:
            # 创建不含等位基因标记的染色体列
            adata_sub.var["chrom_clean"] = adata_sub.var[facet_by].astype(str)
            adata_sub.var["chrom_clean"] = adata_sub.var["chrom_clean"].str.replace("a", "")
            adata_sub.var["chrom_clean"] = adata_sub.var["chrom_clean"].str.replace("b", "")
            
            # 按 (chrom_clean, pos) 分组求均值
            X_data = _get_X(adata_sub, layer)
            if sparse.issparse(X_data):
                X_data = X_data.toarray()
            
            df_data = pd.DataFrame(
                X_data.T,  # 转置：行=变量，列=细胞
                index=adata_sub.var.index,
                columns=adata_sub.obs.index
            )
            
            # 添加分组信息
            df_data["chrom_clean"] = adata_sub.var["chrom_clean"].values
            if "pos" in adata_sub.var.columns:
                df_data["pos"] = adata_sub.var["pos"].values
            else:
                # 从 index 解析位置
                df_data["pos"] = df_data.index.str.split("-").str[1].astype(int)
            
            # 按 (chrom_clean, pos) 分组求均值
            group_cols = ["chrom_clean", "pos"]
            cell_cols = [c for c in df_data.columns if c not in group_cols]
            df_merged = df_data.groupby(group_cols)[cell_cols].mean().reset_index()
            
            # 重建 AnnData
            X_merged = df_merged[cell_cols].values.T  # 转置回来：行=细胞，列=变量
            var_merged = df_merged[["chrom_clean", "pos"]].copy()
            var_merged.index = [f"{row['chrom_clean']}-{row['pos']}" for _, row in var_merged.iterrows()]
            
            # 创建新的 AnnData 对象而不是修改现有的
            import anndata as ad
            adata_plot = ad.AnnData(
                X=X_merged,
                obs=adata_sub.obs.copy(),
                var=var_merged
            )
            adata_plot.var[facet_by] = var_merged["chrom_clean"]
        else:
            adata_plot = adata_sub
        
        # 3. 分组过滤
        if groups is not None:
            group_mask = adata_plot.obs[group_by].isin(groups)
            adata_plot = adata_plot[group_mask, :]
        
        if adata_plot.n_obs == 0:
            raise ValueError("No observations left after filtering")
        
        # 4. 分面过滤
        if facets is not None:
            facet_mask = adata_plot.var[facet_by].isin(facets)
            adata_plot = adata_plot[:, facet_mask]
        
        if adata_plot.n_vars == 0:
            raise ValueError("No variables left after facet filtering")
        
        # 5. 计算分组均值 - 向量化重构版本
        X_plot = _get_X(adata_plot, layer)
        if sparse.issparse(X_plot):
            X_plot = X_plot.toarray()
        
        # 创建包含数据和分组信息的DataFrame
        data_df = pd.DataFrame(X_plot, index=adata_plot.obs.index)
        data_df[group_by] = adata_plot.obs[group_by].values
        
        # 使用groupby进行向量化分组均值计算
        group_means_matrix = data_df.groupby(group_by).apply(
            lambda x: x.drop(columns=[group_by]).apply(lambda col: np.nanmean(col))
        )
        
        # 构建结果DataFrame
        plot_data = []
        var_info = adata_plot.var[[facet_by, "pos"]].reset_index(drop=True)
        
        for group_name in group_means_matrix.index:
            group_means = group_means_matrix.loc[group_name]
            # 过滤掉NaN值
            valid_mask = ~np.isnan(group_means)
            if valid_mask.any():
                valid_indices = np.where(valid_mask)[0]
                group_data = pd.DataFrame({
                    "group": group_name,
                    facet_by: var_info.iloc[valid_indices][facet_by].values,
                    "pos": var_info.iloc[valid_indices]["pos"].values,
                    "mean": group_means.iloc[valid_indices].values
                })
                plot_data.append(group_data)
        
        plot_df = pd.concat(plot_data, ignore_index=True) if plot_data else pd.DataFrame()
        
        if plot_df.empty:
            raise ValueError("No valid data for plotting")
        
        # 6. 颜色设置
        groups_available = plot_df["group"].unique()
        if color_palette is None:
            color_palette = dict(zip(groups_available, sns.color_palette("husl", len(groups_available))))
        
        # 按 color_palette 的 key 顺序排列 groups，确保行顺序正确
        if color_palette is not None:
            # 使用 color_palette 中存在且在数据中也存在的 groups
            groups = [g for g in color_palette.keys() if g in groups_available]
            # 添加数据中存在但 color_palette 中没有的 groups
            missing_groups = [g for g in groups_available if g not in color_palette]
            groups.extend(missing_groups)
            # 为缺失的 groups 添加颜色
            if missing_groups:
                additional_colors = sns.color_palette("husl", len(missing_groups))
                for i, g in enumerate(missing_groups):
                    color_palette[g] = additional_colors[i]
        else:
            groups = list(groups_available)
        
        # 7. 分面设置
        facet_levels = sorted(plot_df[facet_by].unique(), key=natural_chr_key)
        
        # 计算宽度比例（基于每个分面的数据点数量）
        width_ratios = []
        for facet in facet_levels:
            facet_data = plot_df[plot_df[facet_by] == facet]
            if not facet_data.empty:
                width_ratios.append(max(1, len(facet_data) // 100))  # 最小宽度为1
            else:
                width_ratios.append(1)
        
        # 8. 绘图
        if groups_as_rows:
            nrows, ncols = len(groups), len(facet_levels)
            fig, axes = plt.subplots(
                nrows=nrows, ncols=ncols, figsize=figsize, dpi=dpi,
                sharey=True, squeeze=False,
                gridspec_kw={'wspace': 0.0, 'hspace': 0.0, 'width_ratios': width_ratios}
            )
            fig.subplots_adjust(left=0.10, right=0.995, top=0.995, bottom=0.14,
                                wspace=0.0, hspace=0.0)

            for r, g in enumerate(groups):
                col = color_palette[g]
                for c, ch in enumerate(facet_levels):
                    ax = axes[r, c]
                    d = plot_df[(plot_df["group"] == g) & (plot_df[facet_by] == ch)]
                    if not d.empty:
                        ax.set_xlim(d["pos"].min(), d["pos"].max())
                        ax.scatter(d["pos"].values, d["mean"].values, s=s, alpha=alpha, c=[col])
                    
                    # 设置y轴范围
                    if ylim is not None:
                        ax.set_ylim(ylim)

                    ax.grid(True, linestyle="--", alpha=0.25)

                    # X 轴隐藏，用文本写在最下行
                    ax.set_xticks([])
                    ax.tick_params(axis="x", length=0)
                    if r == nrows - 1:
                        ax.text(0.5, -0.16, ch.replace("chr", "Chr "),
                                transform=ax.transAxes, ha="center", va="top",rotation=25)

                    # Y 轴：首列显示，其它列隐藏外观
                    if c == 0:
                        ax.tick_params(axis="y", left=True, labelleft=True, length=2, labelsize=8)
                    else:
                        ax.tick_params(axis="y", left=False, labelleft=False, length=0)

                    # 行标（分组名）
                    if c == 0:
                        ax.text(-0.06, 0.5, str(g), transform=ax.transAxes,
                                va='center', ha='right')

            return fig

        else:
            # 单行版本：不同分组叠加到同一行各分面
            ncols = len(facet_levels)
            fig, axes = plt.subplots(
                nrows=1, ncols=ncols, figsize=figsize, dpi=dpi,
                sharey=True, squeeze=False,
                gridspec_kw={'wspace': 0.0, 'hspace': 0.0, 'width_ratios': width_ratios}
            )
            axrow = axes[0]
            fig.subplots_adjust(left=0.10, right=0.995, top=0.995, bottom=0.14,
                                wspace=0.0, hspace=0.0)

            for c, ch in enumerate(facet_levels):
                ax = axrow[c]
                sub = plot_df[plot_df[facet_by] == ch]
                for g in groups:
                    d = sub[sub["group"] == g]
                    if not d.empty:
                        ax.set_xlim(d["pos"].min(), d["pos"].max())
                        ax.scatter(d["pos"].values, d["mean"].values, s=s, alpha=alpha, c=[color_palette[g]])
                
                # 设置y轴范围
                if ylim is not None:
                    ax.set_ylim(ylim)
                    
                ax.grid(True, linestyle="--", alpha=0.25)

                # x 轴用文本
                ax.set_xticks([])
                ax.tick_params(axis="x", length=0)
                ax.text(0.5, -0.16, ch.replace("chr", "Chr "),
                        transform=ax.transAxes, ha="center", va="top",rotation=25)

                # y 轴：仅第一列显示
                if c == 0:
                    ax.tick_params(axis="y", left=True, labelleft=True, length=2, labelsize=8)
                else:
                    ax.tick_params(axis="y", left=False, labelleft=False, length=0)

            # 图例
            if legend_outside:
                from matplotlib.lines import Line2D
                handles = [Line2D([0], [0], marker='o', linestyle='',
                                  markersize=np.sqrt(s), color=color_palette[g]) for g in groups]
                labels = list(groups)
                fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.01, 0.5), title=group_by)

            return fig
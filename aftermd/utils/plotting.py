import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union


class PlotManager:
    """Plotting utility for visualizing MD analysis results."""
    
    def __init__(self, style: str = "seaborn-v0_8", figsize: Tuple[int, int] = (10, 6)):
        """
        Initialize PlotManager.
        
        Args:
            style: Matplotlib style
            figsize: Default figure size
        """
        plt.style.use(style)
        self.figsize = figsize
        self.colors = sns.color_palette("husl", 10)
        
    def read_xvg_file(self, file_path: str) -> pd.DataFrame:
        """
        Read GROMACS XVG file and return as DataFrame.
        
        Args:
            file_path: Path to XVG file
            
        Returns:
            DataFrame with time series data
        """
        data = []
        headers = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('@') and 'legend' in line:
                    legend = line.split('"')[1]
                    headers.append(legend)
                elif line.startswith('#') or line.startswith('@'):
                    continue
                else:
                    values = [float(x) for x in line.split()]
                    data.append(values)
        
        df = pd.DataFrame(data)
        if headers and len(headers) == len(df.columns) - 1:
            df.columns = ['Time'] + headers
        elif len(df.columns) == 2:
            df.columns = ['Time', 'Value']
        
        return df
    
    def plot_time_series(self, 
                        data: Union[str, pd.DataFrame],
                        title: str = "Time Series",
                        xlabel: str = "Time (ps)",
                        ylabel: str = "Value",
                        output_path: Optional[str] = None,
                        interactive: bool = False) -> None:
        """
        Plot time series data.
        
        Args:
            data: XVG file path or DataFrame
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            output_path: Path to save plot
            interactive: Use plotly for interactive plot
        """
        if isinstance(data, str):
            df = self.read_xvg_file(data)
        else:
            df = data
            
        if interactive:
            self._plot_interactive_time_series(df, title, xlabel, ylabel, output_path)
        else:
            self._plot_static_time_series(df, title, xlabel, ylabel, output_path)
    
    def _plot_static_time_series(self, df: pd.DataFrame, title: str, 
                               xlabel: str, ylabel: str, output_path: Optional[str]):
        """Create static time series plot with matplotlib."""
        fig, ax = plt.subplots(figsize=self.figsize)
        
        time_col = df.columns[0]
        for i, col in enumerate(df.columns[1:]):
            ax.plot(df[time_col], df[col], label=col, color=self.colors[i % len(self.colors)])
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        
        if len(df.columns) > 2:
            ax.legend()
            
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
        plt.show()
    
    def _plot_interactive_time_series(self, df: pd.DataFrame, title: str,
                                    xlabel: str, ylabel: str, output_path: Optional[str]):
        """Create interactive time series plot with plotly."""
        fig = go.Figure()
        
        time_col = df.columns[0]
        for col in df.columns[1:]:
            fig.add_trace(go.Scatter(
                x=df[time_col],
                y=df[col],
                mode='lines',
                name=col
            ))
        
        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            hovermode='x unified'
        )
        
        if output_path:
            if output_path.endswith('.html'):
                fig.write_html(output_path)
            else:
                fig.write_image(output_path)
        
        fig.show()
    
    def plot_distribution(self, 
                         data: Union[List, np.ndarray, pd.Series],
                         title: str = "Distribution",
                         xlabel: str = "Value",
                         bins: int = 50,
                         output_path: Optional[str] = None) -> None:
        """
        Plot distribution histogram.
        
        Args:
            data: Data to plot
            title: Plot title
            xlabel: X-axis label
            bins: Number of histogram bins
            output_path: Path to save plot
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        ax.hist(data, bins=bins, alpha=0.7, density=True, color=self.colors[0])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Density")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
        plt.show()
    
    def plot_heatmap(self, 
                    data: Union[np.ndarray, pd.DataFrame],
                    title: str = "Heatmap",
                    xlabel: str = "X",
                    ylabel: str = "Y",
                    output_path: Optional[str] = None) -> None:
        """
        Plot 2D heatmap.
        
        Args:
            data: 2D data array or DataFrame
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
            output_path: Path to save plot
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        if isinstance(data, pd.DataFrame):
            sns.heatmap(data, ax=ax, cmap='viridis')
        else:
            im = ax.imshow(data, cmap='viridis', aspect='auto')
            plt.colorbar(im, ax=ax)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
        plt.show()
    
    def plot_multiple_xvg(self, 
                         file_paths: List[str],
                         labels: Optional[List[str]] = None,
                         title: str = "Multiple Time Series",
                         output_path: Optional[str] = None) -> None:
        """
        Plot multiple XVG files on the same plot.
        
        Args:
            file_paths: List of XVG file paths
            labels: Optional labels for each file
            title: Plot title
            output_path: Path to save plot
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        if labels is None:
            labels = [Path(fp).stem for fp in file_paths]
        
        for i, (file_path, label) in enumerate(zip(file_paths, labels)):
            df = self.read_xvg_file(file_path)
            time_col = df.columns[0]
            value_col = df.columns[1]
            
            ax.plot(df[time_col], df[value_col], 
                   label=label, color=self.colors[i % len(self.colors)])
        
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("Value")
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
        plt.show()
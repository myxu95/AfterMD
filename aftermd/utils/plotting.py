import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union

# Optional plotly import
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


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
            
        if interactive and PLOTLY_AVAILABLE:
            self._plot_interactive_time_series(df, title, xlabel, ylabel, output_path)
        elif interactive and not PLOTLY_AVAILABLE:
            print("⚠️ Plotly not available, falling back to static plot")
            self._plot_static_time_series(df, title, xlabel, ylabel, output_path)
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
        if not PLOTLY_AVAILABLE:
            raise ImportError("Plotly is required for interactive plots. Install with: pip install plotly")
            
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

    def plot_rmsd(self,
                  rmsd_file: Union[str, List[str]],
                  title: Optional[str] = None,
                  xlabel: str = "Time (ps)",
                  ylabel: str = "RMSD (nm)",
                  output_path: Optional[str] = None,
                  show_stats: bool = True,
                  moving_average: Optional[int] = None,
                  interactive: bool = False) -> None:
        """
        Plot RMSD data from XVG file(s) with enhanced visualization.

        Args:
            rmsd_file: Path to RMSD XVG file or list of paths for comparison
            title: Plot title (auto-generated if None)
            xlabel: X-axis label
            ylabel: Y-axis label
            output_path: Output file path for saving plot
            show_stats: Show statistics (mean, std) on plot
            moving_average: Window size for moving average (optional)
            interactive: Use interactive plotly plot
        """
        if isinstance(rmsd_file, str):
            rmsd_files = [rmsd_file]
        else:
            rmsd_files = rmsd_file

        if interactive and PLOTLY_AVAILABLE:
            self._plot_rmsd_interactive(rmsd_files, title, xlabel, ylabel,
                                      output_path, show_stats, moving_average)
        else:
            self._plot_rmsd_static(rmsd_files, title, xlabel, ylabel,
                                 output_path, show_stats, moving_average)

    def _plot_rmsd_static(self, rmsd_files: List[str], title: Optional[str],
                         xlabel: str, ylabel: str, output_path: Optional[str],
                         show_stats: bool, moving_average: Optional[int]) -> None:
        """Static RMSD plot using matplotlib."""
        fig, ax = plt.subplots(figsize=self.figsize)

        for i, file_path in enumerate(rmsd_files):
            # Read RMSD data
            df = self.read_xvg_file(file_path)
            if df.empty:
                continue

            time_col = df.columns[0]
            rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

            time = df[time_col].values
            rmsd = df[rmsd_col].values

            # Convert time from ps to ns for better readability
            time_ns = time / 1000

            file_name = Path(file_path).stem
            color = self.colors[i % len(self.colors)]

            # Plot main RMSD curve
            ax.plot(time_ns, rmsd, color=color, linewidth=1.5,
                   label=file_name, alpha=0.8)

            # Add moving average if requested
            if moving_average and len(rmsd) >= moving_average:
                rmsd_ma = pd.Series(rmsd).rolling(window=moving_average, center=True).mean()
                ax.plot(time_ns, rmsd_ma, color=color, linewidth=2.5,
                       linestyle='--', alpha=0.9,
                       label=f'{file_name} (MA {moving_average})')

            # Add statistics if requested
            if show_stats:
                mean_rmsd = np.mean(rmsd)
                std_rmsd = np.std(rmsd)

                # Add horizontal line for mean
                ax.axhline(y=mean_rmsd, color=color, linestyle=':',
                          alpha=0.6, linewidth=1)

                # Add text with statistics
                stats_text = f'{file_name}: μ={mean_rmsd:.2f}±{std_rmsd:.2f} nm'
                ax.text(0.02, 0.98 - i*0.05, stats_text,
                       transform=ax.transAxes, fontsize=9,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor=color, alpha=0.2))

        # Customize plot
        ax.set_xlabel(f"{xlabel.replace('(ps)', '(ns)')}")
        ax.set_ylabel(ylabel)

        if title is None:
            if len(rmsd_files) == 1:
                title = f"RMSD Analysis - {Path(rmsd_files[0]).stem}"
            else:
                title = f"RMSD Comparison ({len(rmsd_files)} trajectories)"
        ax.set_title(title, fontsize=14, fontweight='bold')

        # Add grid and legend
        ax.grid(True, alpha=0.3)
        if len(rmsd_files) > 1 or moving_average:
            ax.legend(loc='best', framealpha=0.9)

        # Set y-axis to start from 0
        ax.set_ylim(bottom=0)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"RMSD plot saved to: {output_path}")

        plt.show()

    def _plot_rmsd_interactive(self, rmsd_files: List[str], title: Optional[str],
                              xlabel: str, ylabel: str, output_path: Optional[str],
                              show_stats: bool, moving_average: Optional[int]) -> None:
        """Interactive RMSD plot using plotly."""
        fig = go.Figure()

        for i, file_path in enumerate(rmsd_files):
            # Read RMSD data
            df = self.read_xvg_file(file_path)
            if df.empty:
                continue

            time_col = df.columns[0]
            rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

            time = df[time_col].values
            rmsd = df[rmsd_col].values

            # Convert time from ps to ns
            time_ns = time / 1000

            file_name = Path(file_path).stem
            color = px.colors.qualitative.Set1[i % len(px.colors.qualitative.Set1)]

            # Add main RMSD trace
            fig.add_trace(go.Scatter(
                x=time_ns, y=rmsd,
                mode='lines',
                name=file_name,
                line=dict(color=color, width=2),
                hovertemplate=f'<b>{file_name}</b><br>Time: %{{x:.1f}} ns<br>RMSD: %{{y:.3f}} nm<extra></extra>'
            ))

            # Add moving average if requested
            if moving_average and len(rmsd) >= moving_average:
                rmsd_ma = pd.Series(rmsd).rolling(window=moving_average, center=True).mean()
                fig.add_trace(go.Scatter(
                    x=time_ns, y=rmsd_ma,
                    mode='lines',
                    name=f'{file_name} (MA {moving_average})',
                    line=dict(color=color, width=3, dash='dash'),
                    hovertemplate=f'<b>{file_name} MA</b><br>Time: %{{x:.1f}} ns<br>RMSD: %{{y:.3f}} nm<extra></extra>'
                ))

            # Add mean line if requested
            if show_stats:
                mean_rmsd = np.mean(rmsd)
                fig.add_hline(y=mean_rmsd, line_dash="dot", line_color=color,
                             annotation_text=f'{file_name} mean: {mean_rmsd:.2f} nm',
                             annotation_position="top left")

        # Update layout
        if title is None:
            if len(rmsd_files) == 1:
                title = f"RMSD Analysis - {Path(rmsd_files[0]).stem}"
            else:
                title = f"RMSD Comparison ({len(rmsd_files)} trajectories)"

        fig.update_layout(
            title=title,
            xaxis_title=xlabel.replace('(ps)', '(ns)'),
            yaxis_title=ylabel,
            hovermode='x unified',
            template='plotly_white',
            width=1000,
            height=600,
            yaxis=dict(rangemode='tozero')  # Start y-axis from 0
        )

        if output_path:
            if output_path.endswith('.html'):
                fig.write_html(output_path)
            else:
                fig.write_image(output_path, width=1000, height=600)
            print(f"Interactive RMSD plot saved to: {output_path}")

        fig.show()

    def plot_rmsd_distribution(self,
                              rmsd_file: Union[str, List[str]],
                              title: Optional[str] = None,
                              output_path: Optional[str] = None,
                              bins: int = 50,
                              show_kde: bool = True) -> None:
        """
        Plot RMSD distribution histogram with KDE.

        Args:
            rmsd_file: Path to RMSD XVG file or list of paths
            title: Plot title
            output_path: Output file path
            bins: Number of histogram bins
            show_kde: Show kernel density estimation
        """
        if isinstance(rmsd_file, str):
            rmsd_files = [rmsd_file]
        else:
            rmsd_files = rmsd_file

        fig, ax = plt.subplots(figsize=self.figsize)

        for i, file_path in enumerate(rmsd_files):
            # Read RMSD data
            df = self.read_xvg_file(file_path)
            if df.empty:
                continue

            rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
            rmsd = df[rmsd_col].values

            file_name = Path(file_path).stem
            color = self.colors[i % len(self.colors)]

            # Plot histogram
            ax.hist(rmsd, bins=bins, alpha=0.6, color=color,
                   label=f'{file_name} (n={len(rmsd)})', density=True)

            # Add KDE if requested
            if show_kde:
                try:
                    from scipy.stats import gaussian_kde
                    kde = gaussian_kde(rmsd)
                    x_range = np.linspace(rmsd.min(), rmsd.max(), 200)
                    ax.plot(x_range, kde(x_range), color=color, linewidth=2)
                except ImportError:
                    print("Warning: scipy not available, skipping KDE")

            # Add statistics
            mean_rmsd = np.mean(rmsd)
            ax.axvline(mean_rmsd, color=color, linestyle='--', linewidth=2,
                      label=f'{file_name} mean: {mean_rmsd:.2f} nm')

        # Customize plot
        ax.set_xlabel("RMSD (nm)")
        ax.set_ylabel("Density")

        if title is None:
            if len(rmsd_files) == 1:
                title = f"RMSD Distribution - {Path(rmsd_files[0]).stem}"
            else:
                title = f"RMSD Distribution Comparison"
        ax.set_title(title, fontsize=14, fontweight='bold')

        ax.grid(True, alpha=0.3)
        ax.legend()

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"RMSD distribution plot saved to: {output_path}")

        plt.show()

    def plot_rmsd_convergence(self,
                             rmsd_file: str,
                             window_sizes: List[int] = [100, 500, 1000, 2000],
                             title: Optional[str] = None,
                             output_path: Optional[str] = None) -> None:
        """
        Plot RMSD convergence analysis with different time windows.

        Args:
            rmsd_file: Path to RMSD XVG file
            window_sizes: List of window sizes for moving average
            title: Plot title
            output_path: Output file path
        """
        # Read RMSD data
        df = self.read_xvg_file(rmsd_file)
        if df.empty:
            print(f"Warning: Could not read data from {rmsd_file}")
            return

        time_col = df.columns[0]
        rmsd_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

        time = df[time_col].values / 1000  # Convert to ns
        rmsd = df[rmsd_col].values

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(self.figsize[0], self.figsize[1]*1.5))

        # Plot 1: Raw RMSD with different moving averages
        ax1.plot(time, rmsd, color='lightgray', alpha=0.6, linewidth=0.8, label='Raw RMSD')

        for i, window in enumerate(window_sizes):
            if len(rmsd) >= window:
                rmsd_ma = pd.Series(rmsd).rolling(window=window, center=True).mean()
                color = self.colors[i % len(self.colors)]
                ax1.plot(time, rmsd_ma, color=color, linewidth=2,
                        label=f'MA {window} frames')

        ax1.set_xlabel("Time (ns)")
        ax1.set_ylabel("RMSD (nm)")
        ax1.set_title("RMSD with Moving Averages", fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax1.set_ylim(bottom=0)

        # Plot 2: Cumulative mean convergence
        cumulative_mean = np.cumsum(rmsd) / np.arange(1, len(rmsd) + 1)
        ax2.plot(time, cumulative_mean, color='darkblue', linewidth=2, label='Cumulative Mean')

        # Add horizontal line for final mean
        final_mean = np.mean(rmsd)
        ax2.axhline(final_mean, color='red', linestyle='--', linewidth=2,
                   label=f'Final Mean: {final_mean:.2f} nm')

        # Add convergence threshold (±5% of final mean)
        threshold = 0.05 * final_mean
        ax2.fill_between(time, final_mean - threshold, final_mean + threshold,
                        alpha=0.2, color='red', label='±5% threshold')

        ax2.set_xlabel("Time (ns)")
        ax2.set_ylabel("Cumulative Mean RMSD (nm)")
        ax2.set_title("RMSD Convergence Analysis", fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.set_ylim(bottom=0)

        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        else:
            file_name = Path(rmsd_file).stem
            fig.suptitle(f"RMSD Convergence Analysis - {file_name}",
                        fontsize=16, fontweight='bold')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"RMSD convergence plot saved to: {output_path}")

        plt.show()
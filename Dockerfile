# Use official lightweight Python
FROM python:3.11-slim-bookworm

# Install system packages for Perl, C++ binaries, bioinformatics tools, and Plotly/Kaleido image export
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    wget \
    perl \
    make \
    gcc \
    g++ \
    bcftools \
    samtools \
    tabix \
    chromium \
    && rm -rf /var/lib/apt/lists/*

# Let Kaleido/Plotly find Chromium inside Docker
ENV CHROME_PATH=/usr/bin/chromium
ENV BROWSER_PATH=/usr/bin/chromium

# Tell Python to buffer output
ENV PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Copy requirements first for better Docker cache
COPY requirements.txt .

# Install Python packages
RUN pip install --upgrade "pip<25" wheel \
    && pip install "setuptools==57.5.0" \
    && pip install --no-build-isolation --no-cache-dir -r requirements.txt
# Copy all app files into container
COPY . .

# Make sure Perl scripts & SURVIVOR binary are executable
RUN find EvalSVcallers-master/scripts -name "*.pl" -exec chmod +x {} \; || true
RUN chmod +x SURVIVOR/Debug/SURVIVOR || true

# Create persistent app folders
RUN mkdir -p uploaded_files \
    uploaded_files/survivor_output \
    uploaded_files/evalsvcallers_output \
    uploaded_files/truvari_output \
    uploaded_files/visualization_output \
    uploaded_files/metrics_output \
    uploaded_files/reference_files \
    uploaded_files/example_downloads

# Expose port
EXPOSE 8040

# Run Dash app
CMD ["python", "app.py"]
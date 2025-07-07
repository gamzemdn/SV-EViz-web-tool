# Use official lightweight Python
FROM python:3.11-slim

# Install system packages for Perl and C++ binaries
RUN apt-get update && apt-get install -y \
    perl \
    make \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy all app files into container
COPY . .

# Make sure Perl scripts & SURVIVOR binary are executable
RUN chmod +x EvalSVcallers-master/scripts/*.pl
RUN chmod +x SURVIVOR/Debug/SURVIVOR

# Install Python packages
RUN pip install --no-cache-dir -r requirements.txt

# Expose port
EXPOSE 8040

# Tell Python to buffer output
ENV PYTHONUNBUFFERED=1

# Run your Dash app
CMD ["python", "app.py"]

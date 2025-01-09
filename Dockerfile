FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy application files
COPY . /srv/shiny-server/app
WORKDIR /srv/shiny-server/app

# Install R packages
RUN Rscript requirements.R

# Install Python dependencies
RUN pip3 install -r requirements.txt

# Expose port
EXPOSE 3838

# Run app
CMD ["/usr/bin/shiny-server"]

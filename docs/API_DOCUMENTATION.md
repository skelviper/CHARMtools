# API Documentation

CHARMtools provides a RESTful API built with FastAPI for programmatic access to 3D chromatin analysis capabilities. This API enables integration with web applications, automated workflows, and remote analysis systems.

## Overview

The CHARMtools API provides:
- RESTful endpoints for all major functionality
- Interactive documentation with Swagger UI
- JSON-based request/response format
- Authentication and rate limiting
- Batch processing capabilities

## Getting Started

### Starting the API Server

```bash
# Using the provided startup script
python start_api.py

# Or directly with uvicorn
uvicorn CHARMtools.api:app --host 0.0.0.0 --port 8000

# With custom configuration
python start_api.py --host 127.0.0.1 --port 8080 --workers 4
```

### Accessing Documentation

Once the server is running, access the interactive documentation:

- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc
- **OpenAPI Schema**: http://localhost:8000/openapi.json

## API Endpoints

### Core Endpoints

#### Health Check
```http
GET /health
```

**Response**:
```json
{
  "status": "healthy",
  "version": "1.0.0",
  "timestamp": "2024-01-01T12:00:00Z"
}
```

#### Root Information
```http
GET /
```

**Response**:
```json
{
  "message": "CHARMtools API",
  "version": "1.0.0",
  "documentation": "/docs",
  "endpoints": ["/cell3d", "/analysis", "/visualization"]
}
```

### Cell3D Management

#### Create Cell3D Object
```http
POST /cell3d/create
```

**Request Body**:
```json
{
  "tdg_path": "/path/to/structure.3dg",
  "pairs_path": "/path/to/contacts.pairs.gz",
  "on_disk": false,
  "name": "sample_cell"
}
```

**Response**:
```json
{
  "cell_id": "cell_12345",
  "status": "created",
  "info": {
    "name": "sample_cell",
    "chromosomes": ["chr1", "chr2", "chr3"],
    "total_beads": 50000,
    "memory_usage": "245MB"
  }
}
```

#### Load Existing Cell3D
```http
POST /cell3d/load
```

**Request Body**:
```json
{
  "file_path": "/path/to/saved_cell.h5",
  "format": "hdf5"
}
```

#### Get Cell3D Information
```http
GET /cell3d/{cell_id}/info
```

**Response**:
```json
{
  "cell_id": "cell_12345",
  "name": "sample_cell",
  "chromosomes": ["chr1", "chr2", "chr3"],
  "total_beads": 50000,
  "features": ["enhancers", "promoters"],
  "memory_usage": "245MB",
  "created_at": "2024-01-01T12:00:00Z"
}
```

#### List All Cell3D Objects
```http
GET /cell3d/list
```

**Response**:
```json
{
  "cells": [
    {
      "cell_id": "cell_12345",
      "name": "sample_cell",
      "status": "active",
      "created_at": "2024-01-01T12:00:00Z"
    }
  ],
  "total": 1
}
```

#### Delete Cell3D Object
```http
DELETE /cell3d/{cell_id}
```

**Response**:
```json
{
  "message": "Cell3D object deleted successfully",
  "cell_id": "cell_12345"
}
```

### Feature Management

#### Add Genomic Features
```http
POST /cell3d/{cell_id}/features/add
```

**Request Body**:
```json
{
  "feature_name": "enhancers",
  "bed_file_path": "/path/to/enhancers.bed",
  "feature_type": "genomic_regions"
}
```

**Response**:
```json
{
  "message": "Feature added successfully",
  "feature_name": "enhancers",
  "regions_added": 1250
}
```

### Analysis Endpoints

#### Spatial Analysis
```http
POST /cell3d/{cell_id}/analysis/spatial
```

**Request Body**:
```json
{
  "analysis_type": "distance_analysis",
  "parameters": {
    "features": ["enhancers", "promoters"],
    "distance_threshold": 50000,
    "method": "euclidean"
  }
}
```

**Response**:
```json
{
  "analysis_id": "analysis_67890",
  "status": "completed",
  "results": {
    "mean_distance": 25.7,
    "median_distance": 22.1,
    "contact_probability": 0.15
  }
}
```

#### Structure Analysis
```http
POST /cell3d/{cell_id}/analysis/structure
```

**Request Body**:
```json
{
  "analysis_type": "compactness",
  "chromosomes": ["chr1", "chr2"],
  "parameters": {
    "method": "radius_of_gyration"
  }
}
```

### Visualization Endpoints

#### Generate 3D Plot
```http
POST /cell3d/{cell_id}/visualization/3d
```

**Request Body**:
```json
{
  "plot_type": "structure_3d",
  "parameters": {
    "chromosomes": ["chr1"],
    "color_by": "chromosome",
    "format": "html",
    "interactive": true
  }
}
```

**Response**:
```json
{
  "plot_id": "plot_11111",
  "status": "generated",
  "download_url": "/downloads/plot_11111.html",
  "preview_url": "/preview/plot_11111"
}
```

#### Generate Contact Matrix
```http
POST /cell3d/{cell_id}/visualization/matrix
```

**Request Body**:
```json
{
  "chromosome": "chr1",
  "resolution": 40000,
  "region": "chr1:1000000-5000000",
  "format": "png",
  "colormap": "Reds"
}
```

### Data Export

#### Export Cell3D Data
```http
POST /cell3d/{cell_id}/export
```

**Request Body**:
```json
{
  "format": "hdf5",
  "include_features": true,
  "include_analysis": false,
  "compression": "gzip"
}
```

**Response**:
```json
{
  "export_id": "export_22222",
  "status": "completed",
  "download_url": "/downloads/export_22222.h5",
  "file_size": "125MB"
}
```

## Request/Response Models

### Cell3D Creation Request
```json
{
  "tdg_path": "string (optional)",
  "pairs_path": "string (optional)",
  "on_disk": "boolean (default: false)",
  "name": "string (optional)",
  "metadata": "object (optional)"
}
```

### Feature Addition Request
```json
{
  "feature_name": "string (required)",
  "bed_file_path": "string (optional)",
  "feature_data": "array (optional)",
  "feature_type": "string (default: 'genomic_regions')"
}
```

### Analysis Request
```json
{
  "analysis_type": "string (required)",
  "parameters": "object (required)",
  "chromosomes": "array (optional)",
  "output_format": "string (default: 'json')"
}
```

### Visualization Request
```json
{
  "plot_type": "string (required)",
  "parameters": "object (required)",
  "format": "string (default: 'png')",
  "interactive": "boolean (default: false)"
}
```

## Usage Examples

### Python Client

```python
import requests
import json

# Base URL
base_url = "http://localhost:8000"

# Create Cell3D object
response = requests.post(
    f"{base_url}/cell3d/create",
    json={
        "tdg_path": "/data/structure.3dg",
        "name": "my_cell"
    }
)
cell_data = response.json()
cell_id = cell_data["cell_id"]

# Add features
response = requests.post(
    f"{base_url}/cell3d/{cell_id}/features/add",
    json={
        "feature_name": "enhancers",
        "bed_file_path": "/data/enhancers.bed"
    }
)

# Perform spatial analysis
response = requests.post(
    f"{base_url}/cell3d/{cell_id}/analysis/spatial",
    json={
        "analysis_type": "distance_analysis",
        "parameters": {
            "features": ["enhancers"],
            "distance_threshold": 50000
        }
    }
)
analysis_results = response.json()

# Generate visualization
response = requests.post(
    f"{base_url}/cell3d/{cell_id}/visualization/3d",
    json={
        "plot_type": "structure_3d",
        "parameters": {
            "chromosomes": ["chr1"],
            "color_by": "enhancers",
            "format": "html"
        }
    }
)
plot_info = response.json()
print(f"Plot available at: {plot_info['download_url']}")
```

### cURL Examples

```bash
# Health check
curl -X GET "http://localhost:8000/health"

# Create Cell3D object
curl -X POST "http://localhost:8000/cell3d/create" \
     -H "Content-Type: application/json" \
     -d '{
       "tdg_path": "/data/structure.3dg",
       "name": "test_cell"
     }'

# Get cell information
curl -X GET "http://localhost:8000/cell3d/cell_12345/info"

# Perform analysis
curl -X POST "http://localhost:8000/cell3d/cell_12345/analysis/spatial" \
     -H "Content-Type: application/json" \
     -d '{
       "analysis_type": "distance_analysis",
       "parameters": {
         "distance_threshold": 50000
       }
     }'
```

### JavaScript/Node.js

```javascript
const axios = require('axios');

const baseURL = 'http://localhost:8000';
const client = axios.create({ baseURL });

async function analyzeCell() {
  try {
    // Create Cell3D object
    const createResponse = await client.post('/cell3d/create', {
      tdg_path: '/data/structure.3dg',
      name: 'js_cell'
    });
    
    const cellId = createResponse.data.cell_id;
    console.log(`Created cell: ${cellId}`);
    
    // Add features
    await client.post(`/cell3d/${cellId}/features/add`, {
      feature_name: 'enhancers',
      bed_file_path: '/data/enhancers.bed'
    });
    
    // Perform analysis
    const analysisResponse = await client.post(
      `/cell3d/${cellId}/analysis/spatial`,
      {
        analysis_type: 'distance_analysis',
        parameters: {
          features: ['enhancers'],
          distance_threshold: 50000
        }
      }
    );
    
    console.log('Analysis results:', analysisResponse.data);
    
  } catch (error) {
    console.error('Error:', error.response?.data || error.message);
  }
}

analyzeCell();
```

## Error Handling

### HTTP Status Codes

- **200**: Success
- **201**: Created
- **400**: Bad Request (invalid parameters)
- **404**: Not Found (cell_id not found)
- **422**: Validation Error
- **500**: Internal Server Error

### Error Response Format

```json
{
  "error": {
    "code": "VALIDATION_ERROR",
    "message": "Invalid parameter: distance_threshold must be positive",
    "details": {
      "field": "distance_threshold",
      "value": -100,
      "constraint": "must be > 0"
    }
  }
}
```

### Common Error Codes

- `CELL_NOT_FOUND`: Specified cell_id does not exist
- `VALIDATION_ERROR`: Request validation failed
- `FILE_NOT_FOUND`: Specified file path does not exist
- `ANALYSIS_FAILED`: Analysis computation failed
- `MEMORY_ERROR`: Insufficient memory for operation
- `FORMAT_ERROR`: Invalid file format

## Authentication and Security

### API Key Authentication (Optional)

```python
# Include API key in headers
headers = {
    "X-API-Key": "your-api-key-here",
    "Content-Type": "application/json"
}

response = requests.post(
    f"{base_url}/cell3d/create",
    json=data,
    headers=headers
)
```

### Rate Limiting

- Default: 100 requests per minute per IP
- Configurable via environment variables
- Rate limit headers included in responses

## Performance Considerations

### Async Operations

Long-running operations (analysis, visualization) are handled asynchronously:

```python
# Start analysis
response = requests.post(f"{base_url}/cell3d/{cell_id}/analysis/spatial", json=params)
analysis_id = response.json()["analysis_id"]

# Check status
while True:
    status_response = requests.get(f"{base_url}/analysis/{analysis_id}/status")
    status = status_response.json()["status"]
    
    if status == "completed":
        results = requests.get(f"{base_url}/analysis/{analysis_id}/results")
        break
    elif status == "failed":
        error = requests.get(f"{base_url}/analysis/{analysis_id}/error")
        break
    
    time.sleep(5)  # Wait 5 seconds before checking again
```

### Batch Processing

```python
# Process multiple cells
cells_data = [
    {"tdg_path": "/data/cell1.3dg", "name": "cell1"},
    {"tdg_path": "/data/cell2.3dg", "name": "cell2"},
    {"tdg_path": "/data/cell3.3dg", "name": "cell3"}
]

# Create cells in batch
batch_response = requests.post(
    f"{base_url}/cell3d/batch/create",
    json={"cells": cells_data}
)

cell_ids = [cell["cell_id"] for cell in batch_response.json()["cells"]]
```

## Configuration

### Environment Variables

```bash
# Server configuration
export CHARMTOOLS_HOST=0.0.0.0
export CHARMTOOLS_PORT=8000
export CHARMTOOLS_WORKERS=4

# Security
export CHARMTOOLS_API_KEY=your-secret-key
export CHARMTOOLS_RATE_LIMIT=100

# Storage
export CHARMTOOLS_DATA_DIR=/data
export CHARMTOOLS_TEMP_DIR=/tmp/charmtools
export CHARMTOOLS_MAX_FILE_SIZE=1GB

# Performance
export CHARMTOOLS_MAX_MEMORY=8GB
export CHARMTOOLS_CACHE_SIZE=1GB
```

### Configuration File

```yaml
# config.yaml
server:
  host: 0.0.0.0
  port: 8000
  workers: 4
  reload: false

security:
  api_key: your-secret-key
  rate_limit: 100
  cors_origins: ["*"]

storage:
  data_dir: /data
  temp_dir: /tmp/charmtools
  max_file_size: 1GB

performance:
  max_memory: 8GB
  cache_size: 1GB
  timeout: 300
```

## Deployment

### Docker Deployment

```dockerfile
# Dockerfile
FROM python:3.9-slim

WORKDIR /app
COPY requirements_api.txt .
RUN pip install -r requirements_api.txt

COPY . .
EXPOSE 8000

CMD ["uvicorn", "CHARMtools.api:app", "--host", "0.0.0.0", "--port", "8000"]
```

```bash
# Build and run
docker build -t charmtools-api .
docker run -p 8000:8000 -v /data:/data charmtools-api
```

### Production Deployment

```bash
# Using gunicorn for production
gunicorn CHARMtools.api:app \
  --workers 4 \
  --worker-class uvicorn.workers.UvicornWorker \
  --bind 0.0.0.0:8000 \
  --access-logfile - \
  --error-logfile -
```

## Monitoring and Logging

### Health Monitoring

```python
# Health check endpoint provides system status
response = requests.get(f"{base_url}/health")
health_data = response.json()

print(f"Status: {health_data['status']}")
print(f"Memory usage: {health_data['memory_usage']}")
print(f"Active cells: {health_data['active_cells']}")
```

### Logging Configuration

```python
# Configure logging
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('charmtools_api.log'),
        logging.StreamHandler()
    ]
)
```

## Troubleshooting

### Common Issues

1. **Server won't start**:
   - Check port availability
   - Verify dependencies are installed
   - Check file permissions

2. **Memory errors**:
   - Reduce data size
   - Increase system memory
   - Use on_disk mode for large datasets

3. **File not found errors**:
   - Verify file paths are absolute
   - Check file permissions
   - Ensure files exist on server filesystem

### Debug Mode

```bash
# Start server in debug mode
uvicorn CHARMtools.api:app --reload --log-level debug
```

## API Versioning

The API uses semantic versioning:
- **v1.0.0**: Initial release
- **v1.1.0**: Added batch processing
- **v1.2.0**: Enhanced visualization options

Version information is available at:
```http
GET /version
```

## Support

For API-specific support:
- **Documentation**: Interactive docs at `/docs`
- **GitHub Issues**: Report bugs and feature requests
- **Email**: Contact the development team
- **Community**: Join our discussion forums

---

*This API documentation is automatically generated and kept up-to-date with the latest CHARMtools release. For more information, visit the [CHARMtools documentation](Home.md).*
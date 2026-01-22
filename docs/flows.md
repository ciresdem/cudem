```mermaid
flowchart TD
    %% Define Styles
    classDef external fill:#f9f,stroke:#333,stroke-width:2px;
    classDef module fill:#e1f5fe,stroke:#01579b,stroke-width:2px;
    classDef storage fill:#fff9c4,stroke:#fbc02d,stroke-width:2px;
    classDef product fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px,rx:10,ry:10;
    classDef shared fill:#e0e0e0,stroke:#616161,stroke-width:1px,stroke-dasharray: 5 5;

    %% External Sources
    Sources[("☁️ External Repositories
    NOAA, USGS, NASA, OSM")]:::external

    %% Shared Libraries Subgraph
    subgraph Shared_Libs ["Shared Core Libraries"]
        direction TB
        Regions["Regions.py
        (Bounding Boxes)"]:::shared
        Utils["Utils.py
        (Logging/Progress)"]:::shared
        VDatums["VDatums
        (Vertical Transformation)"]:::shared
    end

    %% Step 1: Acquisition
    subgraph Step1 ["Step 1: Acquisition"]
        Fetches[["Fetches
        (Downloader)"]]:::module
    end

    %% Step 2: Management
    subgraph Step2 ["Step 2: Management"]
        RawFiles[("📂 Local Raw Files
        (Tiff, XYZ, LAS, SHP)")]:::storage
        Dlim[["Dlim
        (Indexer & Processor)"]]:::module
        Datalist[("📄 Datalist File
        (Unified Index)")]:::storage
    end

    %% Step 3: Generation
    subgraph Step3 ["Step 3: Gridding"]
        Waffles[["Waffles
        (Gridding Engine)"]]:::module
    end

    %% Step 4: Product
    subgraph Step4 ["Step 4: Output"]
        DEM([Final DEM / Grid]):::product
        Pancake[["Pancake
        (3D Viz/Conversion)"]]:::module
    end

    %% Connections
    Sources ==> Fetches
    Fetches --> RawFiles
    RawFiles --> Dlim
    Dlim -- "Filters & Datum Shifts" --> Datalist
    Datalist ==> Waffles
    Waffles -- "Interpolation & Clipping" ==> DEM
    DEM -.-> Pancake

    %% Link Shared Libs to Modules (Dotted lines implies dependency)
    Regions -.-> Fetches
    Regions -.-> Dlim
    Regions -.-> Waffles
    VDatums -.-> Dlim
    Utils -.-> Fetches & Dlim & Waffles
    ```
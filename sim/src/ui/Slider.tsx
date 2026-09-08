type Props = {
  label: string;
  value: number;
  min: number;
  max: number;
  step?: number;
  format?: (v: number) => string;
  onChange: (v: number) => void;
};

export function Slider({
  label,
  value,
  min,
  max,
  step = 0.01,
  format = (v) => v.toFixed(3),
  onChange,
}: Props) {
  return (
    <div className="control">
      <label>
        <span>{label}</span>
        <span className="val">{format(value)}</span>
      </label>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(Number(e.target.value))}
      />
    </div>
  );
}
